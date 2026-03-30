import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import pyarrow as pa
import pyarrow.parquet as pq

# 绝对引入根目录库
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.glycan_topology import find_mapped_sugar_units
from lib.cleavage_engine import cleave_glycan_aglycan
from lib.secondary_fragments import extract_secondary_motifs

def process_chunk(chunk):
    """
    处理单一 Chunk 数据流的 6 步内循环逻辑
    """
    # Step 1: 基础清洗 - 过滤 contains_sugar == True
    if 'contains_sugar' in chunk.columns:
        chunk = chunk[chunk['contains_sugar'].astype(str).str.lower() == 'true'].copy()
        
    # Step 1 延续: 依据 standard_inchi_key 去重 (本批次内)
    if 'standard_inchi_key' in chunk.columns:
        chunk = chunk.drop_duplicates(subset=['standard_inchi_key'])
        
    if chunk.empty:
        return chunk

    # 初始化将被注入的新特征列
    chunk['Glycan_SMILES'] = "NULL"
    chunk['Aglycan_SMILES'] = "NULL"
    chunk['NUCLEOTIDES_SMILES'] = "NULL"
    chunk['AminoAcid_SMILES'] = "NULL"
    chunk['Sugar_Sequence'] = "NULL"
    chunk['Sugar_Functional_Group'] = "NULL"

    # 执行行级拓扑解析 (Phase 2: 核心切割)
    for idx, row in chunk.iterrows():
        smiles = str(row.get('canonical_smiles', ''))
        if not smiles or smiles == 'nan':
            continue
            
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
            
        # Step 3: 核心切割 (Phase 2) 优先执行以获得 Aglycan 碎片
        try:
            mapped_units = find_mapped_sugar_units(mol)
            if mapped_units:
                g_tmp, ag_tmp = cleave_glycan_aglycan(mol, mapped_units)
                if g_tmp: chunk.at[idx, 'Glycan_SMILES'] = g_tmp
                if ag_tmp: chunk.at[idx, 'Aglycan_SMILES'] = ag_tmp
            else:
                chunk.at[idx, 'Glycan_SMILES'] = "FALSE_POSITIVE"
                chunk.at[idx, 'Aglycan_SMILES'] = "FALSE_POSITIVE"
        except Exception as e:
            chunk.at[idx, 'Glycan_SMILES'] = "ERROR"
            chunk.at[idx, 'Aglycan_SMILES'] = "ERROR"

    # --- 物理清理逻辑 (Physical Cleanup Logic) ---
    # 如果这一批次清理后空了，直接 continue 进入下一批 (If the chunk is empty after cleaning, skip)
    chunk = chunk[~chunk['Glycan_SMILES'].isin(["FALSE_POSITIVE", "ERROR", "NULL", "nan"])]
    if chunk.empty:
        return chunk

    # --- 清理完毕，才允许幸存的真分子进入后续的 Phase 3 和 5 (Only surviving true molecules proceed) ---
    for idx, row in chunk.iterrows():
        aglycan_smi = chunk.at[idx, 'Aglycan_SMILES']
        glycan_smi = chunk.at[idx, 'Glycan_SMILES']
        
        # Step 2: 动态骨架填补 (User's Rule Fallback)
        if 'murcko_framework' in chunk.columns:
            framework = str(row.get('murcko_framework', ''))
            if framework.strip() in ("", "nan", "NULL", "None") or pd.isna(row.get('murcko_framework')):
                # 检测到 NULL，调用 RDKit 根据计算出的 Aglycan 现场填补 Scaffold
                if aglycan_smi and aglycan_smi not in ("NULL", "ERROR", "FALSE_POSITIVE"):
                    try:
                        ag_mol = Chem.MolFromSmiles(aglycan_smi)
                        if ag_mol:
                            scaffold = MurckoScaffold.GetScaffoldForMol(ag_mol)
                            if scaffold.GetNumAtoms() > 0:
                                chunk.at[idx, 'murcko_framework'] = Chem.MolToSmiles(scaffold)
                    except:
                        pass
                        
        # Step 4: 次级扫描 (Phase 3) 
        if aglycan_smi and aglycan_smi not in ("NULL", "ERROR", "FALSE_POSITIVE"):
            try:
                motifs = extract_secondary_motifs(aglycan_smi)
                if motifs.get('NUCLEOTIDES_SMILES'): chunk.at[idx, 'NUCLEOTIDES_SMILES'] = motifs['NUCLEOTIDES_SMILES']
                if motifs.get('AminoAcid_SMILES'): chunk.at[idx, 'AminoAcid_SMILES'] = motifs['AminoAcid_SMILES']
            except:
                pass

        # Step 5: 序列解码 (Phase 5)
        if glycan_smi and glycan_smi not in ("NULL", "ERROR", "FALSE_POSITIVE"):
            try:
                g_mol = Chem.MolFromSmiles(glycan_smi)
                if g_mol:
                    # 假定旧有 pipeline 提供此接口，保护性调用
                    from lib.monosaccharide_identifier import generate_refined_sequence
                    seq_str, mods_str = generate_refined_sequence(g_mol)
                    if seq_str: chunk.at[idx, 'Sugar_Sequence'] = seq_str
                    if mods_str: chunk.at[idx, 'Sugar_Functional_Group'] = mods_str
            except:
                pass
                
    return chunk


def process_database_in_chunks(input_csv, output_parquet, chunk_size=5000):
    """
    全量数据库增量提取批处理引擎 (OOM-Proof Streaming Engine)
    """
    print(f"🚀 启动抗 OOM 批处理引擎，Chunk Size: {chunk_size}")
    print(f"    来源: {input_csv}")
    print(f"    目标: {output_parquet}")
    
    if not os.path.exists(input_csv):
        print(f"❌ 找不到输入文件: {input_csv} (Please ensure database is located here)")
        return

    # 若目标已存在则清理以防类型叠加崩溃
    if os.path.exists(output_parquet):
        print("⚠️ 检测到旧 Parquet 残留，正在清理覆盖...")
        os.remove(output_parquet)

    # 启动内存友好的 Generator
    chunk_iterator = pd.read_csv(input_csv, chunksize=chunk_size, dtype=str, encoding='utf-8-sig', low_memory=False)
    
    writer = None
    total_processed = 0
    total_valid = 0
    
    for i, chunk in enumerate(chunk_iterator):
        print(f"📦 正在处理第 {i+1} 批次数据... (当前流过内存: {len(chunk)} 行)")
        
        # 将该批 5000 行全数压过我们的 6-Step Pipeline
        processed_chunk = process_chunk(chunk)
        
        if processed_chunk.empty:
            print(f"⏭️ 第 {i+1} 批次无符合条件之糖分子，直接跳过。")
            continue
            
        # 确保所有数据列类型为 string 防止抛出 parquet type 错误
        processed_chunk = processed_chunk.astype(str)
            
        # Step 6: 增量落盘 (Pandas -> PyArrow Table -> Parquet Append)
        table = pa.Table.from_pandas(processed_chunk)
        if writer is None:
            writer = pq.ParquetWriter(output_parquet, table.schema)
            
        writer.write_table(table)
        
        total_processed += len(chunk)
        total_valid += len(processed_chunk)
        print(f"✅ 第 {i+1} 批次处理完毕! 释放内存。此批写入含糖产物: {len(processed_chunk)} 条。")

    if writer:
        writer.close()
        print(f"🎉 全部流式批处理完毕！成功落盘总记录数: {total_valid} / {total_processed}")
        print(f"💾 最终数据输出至: {output_parquet}")
    else:
        print("🕳️ 管线跑完，但未抓取或者写入任何数据。")

if __name__ == "__main__":
    # 配置输入输出根路径映射
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    # 指向用户存放数据的目录
    input_file = os.path.join(base_dir, "data", "Coconut.csv")
    output_file = os.path.join(base_dir, "reports", "Coconut_Processed_Phase6.parquet")
    
    print("---------------------------------------------------------")
    print(" NaturalGlycanLipid 7-Phase Streaming Engine Initialized ")
    print("---------------------------------------------------------")
    
    print("👉 核心批处理器架构搭建完成。")
    print("👉 正式流转处理全量 Coconut 数据库...")
    process_database_in_chunks(input_file, output_file, chunk_size=5000)
