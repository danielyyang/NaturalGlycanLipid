# 项目状态报告 — GlycoNP-Pipeline
# Project Status Report — GlycoNP-Pipeline

> 更新日期 / Updated: 2026-03-31
> 审计版本 / Audit Version: Saponin-Visualization-Complete v2

---

## 一、项目概况 (Project Overview)

**GlycoNP-Pipeline** 是面向糖化学与天然产物化学的自动化计算管线：
- 从 COCONUT 数据库 (~94,000 条含糖化合物) 中自动 **识别、切分和注释糖苷结构**
- 为糖基化天然产物提供 **系统化的骨架分析、序列标注和分类体系**
- 当前最终产出: `GlycoNP_Deep_Enriched_v13_Pruned.csv` + `GlycoNP_Saponin_DB_v13.csv`
- Saponin 皂苷子库的 **40 图表可视化套件 (40-chart Visualization Suite)** 已部署开发完成

---

## 二、代码健康状态 (Code Health Status)

### 2026-03-30 代码连贯性审计 (Coherence Audit)

| 检查项 | 状态 | 说明 |
|--------|------|------|
| Import 链路一致性 | ✅ 已修复 | `cip_exo_engine.py` 现在调用 `lib.virtual_demodify` (含保护门) |
| API 版本统一 | ✅ 已修复 | 审计报告已升级到 `identify_monosaccharide_v10()` |
| 返回值解包顺序 | ✅ 已修复 | `get_split_smiles()` 在 V13 管线中的解包顺序已修正 |
| 路径可移植性 | ✅ 已修复 | V12/V13/Saponin 脚本已消除硬编码 `d:\Glycan_Database` |
| CSV 数据正确性 | ✅ 已修复 | V13 Pruned / Final / Saponin 中 Glycan/Aglycon 列已修正 |
| 文本挖掘统一 | ⚠️ 待优化 | 三套独立的文本挖掘正则待合并到 `lib/` |
| 糖苷键检测位置 | ⚠️ 待优化 | `detectAllGlycosidicBonds()` 仍在 `run_v12_full_pipeline.py` 中 |

### 核心依赖链路 (Canonical Import Chain)

```
scripts/run_v13_full_pipeline.py
  ├─→ lib/glycan_topology.py          [糖环检测 + 切分]
  ├─→ lib/monosaccharide_identifier.py [v10 鉴定 + 序列]
  │     ├─→ lib/cip_exo_engine.py     [CIP+Exo 指纹]
  │     │     └─→ lib/virtual_demodify.py  [虚拟脱修饰] ← 修复后
  │     └─→ lib/glycan_reference_library.py [参考库]
  ├─→ lib/feature_extractor.py        [Murcko + FP]
  └─→ scripts/run_v12_full_pipeline.py [extractSugarFromName, ...]
```

**红线规则**: `lib/` 层 **禁止** 调用 `scripts/` 层（已于 2026-03-30 修复违规项）。

---

## 三、已实现功能 (Implemented Features) ✅

### Phase 1-3: 数据初始化 + 结构切分 + 二级扫描
- [x] COCONUT 全量读取 + InChIKey 去重 (94,242 条)
- [x] 糖环识别 (五/六元含氧饱和环 + 拓扑验证)
- [x] 不饱和环拒止 / 多环聚醚铁律 / 宏环豁免 / 缩醛桥豁免
- [x] 精确糖苷/苷元分离 (碳原子守恒断言)
- [x] 核苷酸识别 (碱基+磷酸共存) / 氨基酸-肽键识别

### Phase 4: 分类学填补
- [x] LOTUS InChIKey Block-1 匹配
- [x] 属名字典推断 (240+ 植物/50+ 真菌/36+ 细菌属)
- [x] GBIF 物种匹配 / API 缓存

### Phase 5: 糖序列 + 苷元骨架 (核心创新)
- [x] 三层匹配引擎 (Tier 1→2→3)
- [x] CIP+Exo 指纹鉴定引擎 v2.0 (walkSugarRing + 加权匹配)
- [x] 虚拟脱修饰 (19 SMIRKS + 3 保护门: MAX_CUT_ATOMS, MIN_REMAINING_HEAVY, 糖环完整性)
- [x] 异头碳手性剥离 (C1 程序化松弛, C2-C5 绝对严格)
- [x] NLP 救援 (名称统计先验 + 7 策略)
- [x] 120+ 参考糖全库 (KEGG Authority)
- [x] Morgan FP (2048-bit, r=2) + Murcko 骨架 + 全碳拓扑骨架

### Phase 6-7: 分类 + 可视化
- [x] SMARTS 骨架分类 (甾体/三萜/黄酮/…) + 糖脂捕获 + Tanimoto ≥0.85
- [x] 三色高亮 v2.1 (🔴糖核心 / 🟡修饰基团 / 🔵苷元)
- [x] HTML/Excel 报告 + 糖苷键 α/β 标注 + Bond_Detail JSON

---

## 四、基准测试 (Benchmark Status) 🟢

| 基准集 | 条目 | 通过率 | 引擎延迟 | 最后验证 |
|--------|------|--------|----------|----------|
| **统一基准集** (benchmark_unified.json) | 226 | **100%** | 13.2ms/mol | 2026-03-27 |
| Tier A 天然产物 | 50 | **100%** | - | 2026-03-25 |
| 1K 分层压力测试 | 900 | 0% 崩溃率 | S=6.4/M=19.2/C=65.4ms | 2026-03-26 |

### 已知引擎局限
- C-糖苷 CIP 失真 (Vitexin D-Man→需 D-Glc)
- D/L 手性反转 (部分 PubChem SMILES CIP 编码歧义)
- 季碳 SubstructMatch 局限 (L-Cladinose)
- 核苷 C1 手性无法由模板匹配 (退回 Pen/PenN)

---

## 五、数据产出 (Data Outputs)

### 主数据文件

| 文件 | 大小 | 说明 | 状态 |
|------|------|------|------|
| `GlycoNP_Deep_Enriched_v13_Final.csv` | 228MB | V13 完整版 | ✅ Glycan/Aglycon 已修正 |
| `GlycoNP_Deep_Enriched_v13_Pruned.csv` | 163MB | V13 修剪版 | ✅ Glycan/Aglycon 已修正 |
| `GlycoNP_Saponin_DB_v13.csv` | 52MB | Saponin 子库 | ✅ 已同步修正 |

### 关键列说明

| 列名 | 说明 |
|------|------|
| `Glycan_SMILES` | 糖基部分 SMILES (2026-03-30 修正后为正确值) |
| `Aglycon_SMILES` | 苷元部分 SMILES (2026-03-30 修正后为正确值) |
| `Consensus_Sugar_Sequence` | 最终糖序列 (NLP 验证后) |
| `Glycan_Modifications` | 修饰基团标注 |
| `Aglycone_Linkage_Type` | 糖苷键类型 (α/β-O/N/S) |
| `Glycan-Aglycone_Bond_Detail` | 糖苷键 JSON 详情 |
| `Murcko_Scaffold` | 苷元 Bemis-Murcko 骨架 |
| `Organism_Type` | 生物来源 (Plant/Fungi/Bacteria/Marine/Unknown) |
| `Super_Scaffold_Class` | 天然产物超类 |

---

## 六、开发约定 (Development Conventions)

### 文件组织
- **核心算法**: 放在 `lib/` (禁止内联到 scripts)
- **管线执行**: 放在 `scripts/`
- **一次性诊断**: 放在 `scripts/` 但使用后归档到 `scripts/archive/`
- **测试**: 放在 `tests/`, 使用 `pytest`

### 鉴定引擎调用规范
```python
# ✅ 正确: 使用最终版 API
from lib.monosaccharide_identifier import generate_refined_sequence
seq, mods = generate_refined_sequence(mol)

# ✅ 正确: 便捷包装 (SMILES 输入)
from lib.monosaccharide_identifier import analyze_glycan
seq, mods = analyze_glycan(smiles)

# ❌ 禁止: 直接调用基础版 (仅供 v10 内部使用)
# from lib.monosaccharide_identifier import identify_monosaccharide_v2
```

### Import 路径规范
```python
# ✅ 使用相对路径
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# ❌ 禁止硬编码绝对路径
# sys.path.insert(0, r"d:\Glycan_Database")
```

---

## 七、待解决问题 (Open Issues)

| # | 问题 | 优先级 | 说明 |
|:-:|------|:------:|------|
| 1 | 文本挖掘三套字典未统一 | 🟡 P1 | `run_v12`, `rescue_generic_sugars`, `stereochemistry_rescue` 各有独立实现 |
| 2 | `detectAllGlycosidicBonds()` 未抽到 lib/ | 🟡 P1 | 4 个脚本通过 `from scripts.run_v12_full_pipeline` 调用 |
| 3 | API 连通性 | 🟡 | COCONUT ❌ 404, Wikidata ❌ 403, PubChem ✅, GBIF ✅ |
| 4 | ChEMBL 立体化学救援被禁用 | 🟡 | 需 InChIKey 前缀索引 (30GB SQLite) |
| 5 | ~49.4% 化合物无 organism | 🟡 | API 限制导致 |
| 6 | 过渡脚本归档 | 🟢 P2 | 20+ 个旧版脚本需移到 `scripts/archive/` |

---

## 八、运行验证 (Verification)

```bash
# 1. 运行基准测试
python scripts/run_benchmark_200_tier_a.py

# 2. 运行审计修复验证 + CSV 更新
python scripts/run_audit_fixes.py

# 3. 运行单元测试
python -m pytest tests/ -v --tb=short
```
