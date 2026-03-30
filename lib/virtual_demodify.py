"""
虚拟脱修饰概念验证 (Virtual Demodification Proof of Concept)
============================================================

化学红线 (Chemical Red Line):
  切除修饰基团时, 糖环上的连接氧原子 (Anchor Oxygen) 必须保留为 -OH。
  如果连着氧一起切了, 糖就变成脱氧糖, 这是绝对不允许的。
  
  When cleaving modification groups, the anchor oxygen on the sugar ring
  MUST be preserved as -OH. Cutting the oxygen turns the sugar into a
  deoxy sugar — this is absolutely forbidden.

算法 (Algorithm):
  1. 用 SMARTS 定位修饰基团 (锚点 = 连接在糖环上的 O 或 N)
  2. 从匹配中收集 "修饰侧" 原子 (不包括锚点)
  3. 用 RWMol 删除修饰侧原子, 保留锚点
  4. 锚点补氢 → 变成 -OH 或 -NH2
  5. Sanitize + AssignStereochemistry (重建手性)

Usage:
  python scripts/poc_virtual_demodify.py
"""
import sys
import os
from typing import List, Tuple, Optional, Dict

from rdkit import Chem
from rdkit.Chem import AllChem


# =====================================================================
# SMARTS 定义: 修饰基团模式
# [anchor:1] 是连接在糖环上的 O/N — 这个原子必须保留
# 后面的原子是"修饰侧" — 这些原子需要切除
# =====================================================================

# 设计意图 (Design Intent):
#   每个 SMARTS 的第一个原子 ([O:1] 或 [N:1]) 是"锚点" (anchor)。
#   它直接连在糖环碳上, 是糖本固有结构的一部分。
#   切除修饰后, 锚点变成 -OH (氧) 或 -NH2 (氮)。
#   这样糖环保持化学完整 — 不会变成脱氧糖。

DEMOD_PATTERNS = {
    # === O-修饰 (O-linked modifications) ===
    "O-Sulfate": {
        "smarts": "[O:1]S(=O)(=O)[O;H1,O-]",
        "anchor_type": "O",
        "description": "硫酸酯化 — 切除 -SO3H, 保留 -OH",
    },
    "O-Acetyl": {
        "smarts": "[O:1]C(=O)[CH3]",
        "anchor_type": "O",
        "description": "乙酰化 — 切除 -C(=O)CH3, 保留 -OH",
    },
    "O-Methyl": {
        "smarts": "[O:1][CH3;!$([CH3]C=O)]",
        "anchor_type": "O",
        "description": "甲醚化 — 切除 -CH3, 保留 -OH",
    },
    "O-Formyl": {
        "smarts": "[O:1][CH1](=O)",
        "anchor_type": "O",
        "description": "甲酰化 — 切除 -CHO, 保留 -OH",
    },
    "O-Benzoyl": {
        "smarts": "[O:1]C(=O)c1ccccc1",
        "anchor_type": "O",
        "description": "苯甲酰化 — 切除 -C(=O)Ph, 保留 -OH",
    },
    "O-Galloyl": {
        "smarts": "[O:1]C(=O)c1cc(O)c(O)c(O)c1",
        "anchor_type": "O",
        "description": "没食子酰化 — 切除 galloyl, 保留 -OH",
    },
    "O-Phosphate": {
        "smarts": "[O:1]P(=O)([O;H1,O-])[O;H1,O-]",
        "anchor_type": "O",
        "description": "磷酸酯化 — 切除 -PO3H2, 保留 -OH",
    },
    # === N-修饰 (N-linked modifications) ===
    "N-Sulfate": {
        "smarts": "[N:1]S(=O)(=O)[O;H1,O-]",
        "anchor_type": "N",
        "description": "N-硫酸酯化 — 切除 -SO3H, 保留 -NH",
    },
    "N-Acetyl": {
        "smarts": "[N:1]C(=O)[CH3]",
        "anchor_type": "N",
        "description": "N-乙酰化 — 切除 -C(=O)CH3, 保留 -NH",
        # 注意: 不要为 GlcNAc 类剥离 NAc — 这是识别特征
        # Note: Do NOT strip NAc from GlcNAc — it's an identity feature
    },
    "N-Dimethyl": {
        "smarts": "[N:1]([CH3])[CH3]",
        "anchor_type": "N",
        "description": "N-二甲基 — 切除 2×CH3, 保留 -NH (Desosamine 等)",
        # Design Intent: Desosamine has -N(CH3)2 at C3. Strip both methyls
        # to reveal bare -NH2 for CIP matching.
    },
}

# 预编译 SMARTS (Precompile SMARTS patterns)
COMPILED_PATTERNS: Dict[str, Tuple] = {}
for modName, modDef in DEMOD_PATTERNS.items():
    pat = Chem.MolFromSmarts(modDef["smarts"])
    if pat is None:
        print(f"  [ERROR] Invalid SMARTS for {modName}: {modDef['smarts']}")
    else:
        COMPILED_PATTERNS[modName] = (pat, modDef["anchor_type"])


# 保护门常量 (Protection Gate Constants)
# 设计意图: 防止虚拟剥离过度切割苷元 (如 Sennoside A 的蒽醌)
# Prevents over-cutting aglycon groups (e.g., Sennoside A's anthraquinone)
MAX_CUT_ATOMS = 10  # 单次切除不超过 10 个重原子 (标准 Bz=7C, Galloyl=10)
MIN_REMAINING_HEAVY = 5  # 切除后至少保留 5 个重原子 (最小产糖骨架)


def virtualDemodify(
    mol: Chem.Mol,
    modNames: Optional[List[str]] = None,
    maxIterPerMod: int = 10,
) -> Tuple[Chem.Mol, List[str]]:
    """虚拟脱修饰: 切除修饰基团, 保留锚点 O/N 为 -OH/-NH2。
    Virtual demodification: remove modification groups, preserve anchor O/N.

    化学红线 (Chemical Red Line):
      - 切除后锚点 O 必须变成 -OH (不能变成脱氧碳)
      - 切除后锚点 N 必须变成 -NH2 (保留氨基特征)
      - 手性中心 (C2-C5) 必须在操作后重新计算

    保护门 (Protection Gate):
      - 单次切除 > MAX_CUT_ATOMS → 中止 (可能是苷元)
      - 切除后无完整糖环 → 回滚 (切碎了糖环)
      - 切除后重原子 < MIN_REMAINING_HEAVY → 中止

    Args:
        mol: 输入分子 (RDKit Mol)
        modNames: 要剥离的修饰列表 (None=全部 O-修饰)
        maxIterPerMod: 每种修饰最多迭代次数 (防止无限循环)

    Returns:
        (cleanMol, detectedMods): 干净分子 + 检测到的修饰列表
    """
    if modNames is None:
        # 默认只剥离 O-修饰 (Default: O-modifications only)
        # 不剥离 N-Acetyl — 它是 GlcNAc/GalNAc 的识别特征
        modNames = [name for name in COMPILED_PATTERNS
                    if name.startswith("O-")]

    detectedMods = []
    currentMol = Chem.RWMol(Chem.Mol(mol))

    for modName in modNames:
        if modName not in COMPILED_PATTERNS:
            continue
        pattern, anchorType = COMPILED_PATTERNS[modName]

        for iteration in range(maxIterPerMod):
            # 每次迭代重新搜索 (索引可能因前次删除而变化)
            # Re-search each iteration (indices shift after deletion)
            matches = currentMol.GetSubstructMatches(pattern)
            if not matches:
                break

            # 取第一个匹配 (Take first match)
            match = matches[0]
            anchorIdx = match[0]  # [O:1] 或 [N:1] — 锚点, 必须保留

            # 收集修饰侧原子 (非锚点的匹配原子)
            # Collect modification-side atoms (all matched atoms except anchor)
            modSideAtoms = set(match[1:])

            # BFS 从修饰侧向外扩展, 不越过锚点
            # BFS expand outward from mod-side, DO NOT cross anchor
            visited = set(modSideAtoms)
            queue = list(modSideAtoms)
            while queue:
                curr = queue.pop(0)
                atom = currentMol.GetAtomWithIdx(curr)
                for nbr in atom.GetNeighbors():
                    nIdx = nbr.GetIdx()
                    if nIdx not in visited and nIdx != anchorIdx:
                        visited.add(nIdx)
                        queue.append(nIdx)

            allModAtoms = visited

            # === 保护门 1: 原子阈值检查 ===
            # Protection Gate 1: Atom count threshold
            # 如果 BFS 扩展收集了太多原子, 说明可能切到了苷元
            if len(allModAtoms) > MAX_CUT_ATOMS:
                break  # 跳过这种修饰的剩余匹配

            # === 保护门 2: 剩余原子检查 ===
            # Protection Gate 2: Remaining atoms check
            remainingHeavy = currentMol.GetNumHeavyAtoms() - len(allModAtoms)
            if remainingHeavy < MIN_REMAINING_HEAVY:
                break

            # === 保护门 3: 预检环完整性 ===
            # Protection Gate 3: Pre-check ring integrity
            # 确保切除后至少还有一个含 O 的 5/6 元环
            remainingAtoms = set(range(currentMol.GetNumAtoms())) - allModAtoms
            ri = currentMol.GetRingInfo()
            hasValidRingAfterCut = False
            for ring in ri.AtomRings():
                if len(ring) in (5, 6):
                    ringSet = set(ring)
                    # 环原子全部在剩余集合中
                    if ringSet.issubset(remainingAtoms):
                        # 且环中恰好有 1 个 O
                        ringO = sum(1 for i in ring
                                    if currentMol.GetAtomWithIdx(i).GetAtomicNum() == 8)
                        if ringO == 1:
                            hasValidRingAfterCut = True
                            break
            if not hasValidRingAfterCut:
                break  # 切除会破坏糖环 — 中止

            # === 关键步骤: 删除修饰原子, 保留锚点 ===
            # Critical: Remove mod atoms, PRESERVE anchor

            # 1. 先断开锚点与修饰侧的所有键
            #    First disconnect anchor from mod-side
            bondsToBreak = []
            anchorAtom = currentMol.GetAtomWithIdx(anchorIdx)
            for nbr in anchorAtom.GetNeighbors():
                if nbr.GetIdx() in allModAtoms:
                    bond = currentMol.GetBondBetweenAtoms(anchorIdx, nbr.GetIdx())
                    if bond:
                        bondsToBreak.append((anchorIdx, nbr.GetIdx()))

            for a1, a2 in bondsToBreak:
                currentMol.RemoveBond(a1, a2)

            # 2. 从大索引到小索引删除修饰原子 (避免索引偏移)
            #    Remove mod atoms from largest to smallest index
            for atomIdx in sorted(allModAtoms, reverse=True):
                currentMol.RemoveAtom(atomIdx)

            # 3. 锚点可能因删除而索引偏移 — 重新定位
            #    Anchor index may shift due to removal — recalculate
            #    (删除索引 > anchorIdx 的原子不影响 anchorIdx)
            #    (删除索引 < anchorIdx 的原子会让 anchorIdx 减 1)
            newAnchorIdx = anchorIdx
            for removedIdx in sorted(allModAtoms):
                if removedIdx < anchorIdx:
                    newAnchorIdx -= 1

            # 4. 补氢: 锚点变成 -OH 或 -NH2
            #    Cap: anchor becomes -OH or -NH2
            try:
                anchorAtomNew = currentMol.GetAtomWithIdx(newAnchorIdx)
                anchorAtomNew.SetNoImplicit(False)
                anchorAtomNew.SetNumExplicitHs(0)
                anchorAtomNew.SetFormalCharge(0)
            except Exception:
                pass  # 安全失败 (Safe fallback)

            detectedMods.append(modName)

    # === 最终清理: Sanitize + 手性重建 ===
    # Final cleanup: Sanitize + stereochemistry recalculation
    try:
        currentMol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(currentMol)
        cleanMol = currentMol.GetMol()
        cleanMol = Chem.RemoveHs(cleanMol)
        Chem.AssignStereochemistry(cleanMol, force=True, cleanIt=True)
    except Exception as e:
        print(f"  [WARN] Sanitize failed: {e}")
        cleanMol = currentMol.GetMol()

    return cleanMol, detectedMods


# =====================================================================
# 测试用例 (Test Cases)
# =====================================================================

def main():
    print("=" * 70)
    print("  Virtual Demodification PoC")
    print("  虚拟脱修饰概念验证")
    print("=" * 70)

    testCases = [
        # (名称, SMILES, 描述, 预期修饰)
        (
            "D-GlcNAc 6-Sulfate (Benchmark #5)",
            "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",
            "N-Acetyl + O-6-Sulfate; 应产出 D-GlcN 裸糖骨架",
            ["O-Sulfate"],
        ),
        (
            "DiAcetyl D-GlcA (Benchmark #1)",
            "O[C@@H]1[C@H](OC(=O)C)[C@@H](OC(=O)C)[C@H](O)[C@@H](C(=O)O)O1",
            "2x O-Acetyl; 应产出 D-GlcA 裸糖骨架",
            ["O-Acetyl", "O-Acetyl"],
        ),
        (
            "Formylated L-Rha (Benchmark #3)",
            "C[C@@H]1OC(O)[C@H](O)[C@H](OC=O)[C@H]1O",
            "1x O-Formyl; 应产出 L-Rha 裸糖骨架",
            ["O-Formyl"],
        ),
        (
            "Fondaparinux 片段: 硫酸化 GlcN",
            "O[C@H]1[C@@H](O)[C@@H](NS(=O)(=O)O)[C@H](O)O[C@@H]1COS(=O)(=O)O",
            "O-6-Sulfate + N-Sulfate 同时存在; O-Sulfate 应剥离, N-Sulfate 可选",
            ["O-Sulfate"],
        ),
        (
            "Methylated L-Ara (Benchmark #4)",
            "O[C@H]1[C@@H](OC)[C@@H](O)[C@@H](O)CO1",
            "1x O-Methyl; 应产出 L-Ara 裸糖骨架",
            ["O-Methyl"],
        ),
    ]

    for name, smi, desc, expectedMods in testCases:
        print(f"\n--- {name} ---")
        print(f"  描述: {desc}")
        print(f"  原始 SMILES: {smi}")

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print("  [ERROR] RDKit 解析失败!")
            continue

        # 计算脱修饰前的原子计数 (Pre-demod atom counts)
        preO = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
        preC = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        preN = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        preTotal = mol.GetNumHeavyAtoms()
        print(f"  原子计数 (Before): C={preC}, N={preN}, O={preO}, Total={preTotal}")

        # 执行虚拟脱修饰 (Run virtual demodification)
        cleanMol, detectedMods = virtualDemodify(mol)

        cleanSmi = Chem.MolToSmiles(cleanMol)
        postO = sum(1 for a in cleanMol.GetAtoms() if a.GetAtomicNum() == 8)
        postC = sum(1 for a in cleanMol.GetAtoms() if a.GetAtomicNum() == 6)
        postN = sum(1 for a in cleanMol.GetAtoms() if a.GetAtomicNum() == 7)
        postTotal = cleanMol.GetNumHeavyAtoms()

        print(f"  检测到修饰: {detectedMods}")
        print(f"  裸糖 SMILES: {cleanSmi}")
        print(f"  原子计数 (After):  C={postC}, N={postN}, O={postO}, Total={postTotal}")

        # === 化学验证 (Chemical validation) ===
        # 1. 锚点 O 保留检查: O 数不应减少超过 (每个 Sulfate 1 个 S + 3 个 O-远端)
        #    实际: 每个 O-Sulfate 切除 S + 3O (远端), 保留 1O (锚点) → O 应减少 3 per sulfate
        #    每个 O-Acetyl 切除 C-C=O-CH3 (3 atoms), 保留 1O (锚点) → O 应减少 1 per acetyl
        print(f"  O 变化: {preO} -> {postO} (diff={preO - postO})")
        print(f"  C 变化: {preC} -> {postC} (diff={preC - postC})")

        # 2. 手性保留检查 (Chirality preservation check)
        Chem.AssignStereochemistry(cleanMol, force=True, cleanIt=True)
        chiralCenters = Chem.FindMolChiralCenters(cleanMol)
        print(f"  手性中心: {chiralCenters}")
        print()

    # === 附加测试: Fondaparinux 全分子 (含 N-Sulfate 剥离) ===
    print("\n" + "=" * 70)
    print("  Extended Test: Fondaparinux 片段 — 全修饰剥离 (O + N)")
    print("=" * 70)

    fondaFragment = "O[C@H]1[C@@H](O)[C@@H](NS(=O)(=O)O)[C@H](O)O[C@@H]1COS(=O)(=O)O"
    mol = Chem.MolFromSmiles(fondaFragment)
    print(f"  原始: {fondaFragment}")

    # 先只剥 O-Sulfate
    cleanMol1, mods1 = virtualDemodify(mol, modNames=["O-Sulfate"])
    print(f"  仅 O-Sulfate 剥离: {Chem.MolToSmiles(cleanMol1)}")
    print(f"  检测修饰: {mods1}")

    # 再剥 O-Sulfate + N-Sulfate
    cleanMol2, mods2 = virtualDemodify(mol, modNames=["O-Sulfate", "N-Sulfate"])
    print(f"  O+N Sulfate 全剥离: {Chem.MolToSmiles(cleanMol2)}")
    print(f"  检测修饰: {mods2}")

    # 手性保留
    Chem.AssignStereochemistry(cleanMol2, force=True, cleanIt=True)
    chiralFinal = Chem.FindMolChiralCenters(cleanMol2)
    print(f"  手性中心: {chiralFinal}")

    # 裸糖匹配测试 (Bare sugar matching test)
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    try:
        from lib.monosaccharide_identifier import identify_monosaccharide_v2
        from lib.glycan_topology import find_mapped_sugar_units
        units = find_mapped_sugar_units(cleanMol2)
        if units:
            ringAtoms = units[0].get("ring_atoms", [])
            name, anomer = identify_monosaccharide_v2(cleanMol2, ringAtoms)
            print(f"  裸糖匹配结果: {name} (anomer={anomer})")
        else:
            print(f"  裸糖匹配: 未找到糖环")
    except Exception as e:
        print(f"  裸糖匹配失败: {e}")

    print(f"\n{'='*70}")
    print("  PoC Complete")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
