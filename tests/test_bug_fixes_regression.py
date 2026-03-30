"""
回归测试: 修复 5 个管线 BUG / Regression Tests for 5 Pipeline Bug Fixes
[TEST DATA ONLY]

测试用例 (Test Cases):
  1. 硫酸基应归入糖区 (Sulfate atoms should be in glycan zone)
  2. 有苷元的化合物 cleaveAndCap 应返回非空苷元 (Aglycon should not be None)
  3. 超长碳链 (>10 原子) 应归入苷元 (Long chains → aglycon)
  4. D-Rib 环不应被分类为 Non (D-Rib should not be classified as Non)
  5. _countFragmentCarbons 应有深度限制 (BFS should have depth limit)
"""
import os
import sys
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

from rdkit import Chem


# =====================================================================
# Bug 1: 硫酸基/磷酸基应归入糖区 (Sulfate/Phosphate in glycan zone)
# =====================================================================

class TestSulfateColoring:
    """修饰基团 (Sulfate, Phosphate) 应在糖区, 可视化时染黄色而非蓝色。
    Modification groups should be in glycan zone, colored yellow not blue."""

    def test_sulfate_glucosamine_in_glycan_zone(self):
        """6-O-Sulfate-GlcNAc: 所有 S 和 =O 原子应在 glycanAtoms 中"""
        from molecular_visualizer import identifySugarAtomZones
        from glycan_topology import find_mapped_sugar_units

        # 6-O-Sulfate-GlcNAc [TEST DATA ONLY]
        smiles = "O=S(=O)(O)OC[C@H]1OC(O)[C@H](NC(=O)C)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

        units = find_mapped_sugar_units(mol)
        assert len(units) > 0, "Should detect at least 1 sugar unit"

        sugarRingAtoms, substituentAtoms, aglyconAtoms = identifySugarAtomZones(mol, units)
        glycanAtoms = sugarRingAtoms | substituentAtoms

        # 找到所有 S 原子 (Find all S atoms)
        sulfurAtoms = [
            a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 16
        ]
        assert len(sulfurAtoms) > 0, "Should have S atoms"

        # 所有 S 原子应在 glycanAtoms (All S atoms should be in glycan zone)
        for sIdx in sulfurAtoms:
            assert sIdx in glycanAtoms, (
                f"S atom {sIdx} should be in glycanAtoms, not aglyconAtoms"
            )

        # S 原子连接的 =O 也应在 glycanAtoms (S's bonded O should be in glycan)
        for sIdx in sulfurAtoms:
            sAtom = mol.GetAtomWithIdx(sIdx)
            for nbr in sAtom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    assert nbr.GetIdx() in glycanAtoms, (
                        f"O atom {nbr.GetIdx()} bonded to S should be in glycanAtoms"
                    )

    def test_phosphate_mannose_in_glycan_zone(self):
        """Man-6-P: P 和 =O 应在 glycanAtoms 中"""
        from molecular_visualizer import identifySugarAtomZones
        from glycan_topology import find_mapped_sugar_units

        # Man-6-phosphate [TEST DATA ONLY]
        smiles = "O=P(O)(O)OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

        units = find_mapped_sugar_units(mol)
        sugarRingAtoms, substituentAtoms, aglyconAtoms = identifySugarAtomZones(mol, units)
        glycanAtoms = sugarRingAtoms | substituentAtoms

        # P 原子应在 glycanAtoms
        pAtoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
        for pIdx in pAtoms:
            assert pIdx in glycanAtoms, (
                f"P atom {pIdx} should be in glycanAtoms"
            )


# =====================================================================
# Bug 2: cleaveAndCap 应正确返回苷元 (Aglycon should not be None)
# =====================================================================

class TestAglyconPresence:
    """有苷元的化合物 cleaveAndCap 应返回非空结果。
    Compounds with aglycon should return non-None from cleaveAndCap."""

    def test_rutin_has_aglycon(self):
        """芦丁 (Rutin) 应有查耳酮苷元"""
        from molecular_visualizer import identifySugarAtomZones, cleaveAndCap
        from glycan_topology import find_mapped_sugar_units

        RUTIN = "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"
        mol = Chem.MolFromSmiles(RUTIN)
        units = find_mapped_sugar_units(mol)
        sugarRingAtoms, substituentAtoms, aglyconAtoms = identifySugarAtomZones(mol, units)
        glycanAtoms = sugarRingAtoms | substituentAtoms

        aglyconSmi, glycanSmi = cleaveAndCap(
            mol, glycanAtoms, aglyconAtoms,
            minAglyconHeavyAtoms=10, sugarUnits=units,
        )
        assert aglyconSmi is not None, "Rutin aglycon should not be None"
        assert glycanSmi is not None, "Rutin glycan should not be None"

    def test_ginsenoside_has_aglycon(self):
        """人参皂苷 Rg1 应有达玛烷苷元"""
        from molecular_visualizer import identifySugarAtomZones, cleaveAndCap
        from glycan_topology import find_mapped_sugar_units

        RG1 = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(O)C[C@@H](C)O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)C[C@H](O)[C@@H]5[C@@]3(C)CC[C@H](C5(C)C)O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O"
        mol = Chem.MolFromSmiles(RG1)
        units = find_mapped_sugar_units(mol)
        sugarRingAtoms, substituentAtoms, aglyconAtoms = identifySugarAtomZones(mol, units)
        glycanAtoms = sugarRingAtoms | substituentAtoms

        aglyconSmi, glycanSmi = cleaveAndCap(
            mol, glycanAtoms, aglyconAtoms,
            minAglyconHeavyAtoms=10, sugarUnits=units,
        )
        assert aglyconSmi is not None, "Ginsenoside aglycon should not be None"


# =====================================================================
# Bug 3: 修饰基团/苷元阈值 (Modification vs Aglycone Threshold)
# =====================================================================

class TestSubstituentThreshold:
    """超过 10 原子的分支应归入苷元, ≤10 归入修饰。
    Branches >10 atoms → aglycone; ≤10 → substituent."""

    def test_large_branch_is_aglycone(self):
        """自创含 12 原子碳链的糖苷, 该链应归入 aglycone"""
        from glycan_topology import classify_sugar_parts

        # 简单 Glc + 十二碳链 [TEST DATA ONLY]
        # OC[C@H]1OC(OCCCCCCCCCCCC)C(O)C(O)C1O — 糖苷 + C12 链
        smiles = "OC[C@H]1OC(OCCCCCCCCCCCC)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pytest.skip("SMILES parse failed")

        classification = classify_sugar_parts(mol)
        aglyconAtoms = classification['aglycone_atoms']

        # C12 链的原子数 > 10 → 应在 aglycone (12 C + some O)
        assert len(aglyconAtoms) > 0, (
            "Branch with 12 carbons should be classified as aglycone, not substituent"
        )

    def test_small_branch_is_substituent(self):
        """乙酰基 (4 原子) 应归入修饰基团"""
        from glycan_topology import classify_sugar_parts

        # 2,3-di-O-Acetyl-Glucose [TEST DATA ONLY]
        smiles = "OC[C@H]1OC(O)[C@H](OC(=O)C)[C@@H](OC(=O)C)[C@@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pytest.skip("SMILES parse failed")

        classification = classify_sugar_parts(mol)
        substituentAtoms = classification['sugar_substituent_atoms']

        # 乙酰基 (C2H3O ≈ 4 heavy atoms) 应在 substituent
        assert len(substituentAtoms) > 0, (
            "Acetyl groups should be classified as substituents, not aglycone"
        )


# =====================================================================
# Bug 4: D-Rib 不应被分类为 Non (_countFragmentCarbons depth limit)
# =====================================================================

class TestCarbonCountDepthLimit:
    """_countFragmentCarbons 应有深度限制, 防止走入苷元碳骨架。
    _countFragmentCarbons should have depth limit to prevent aglycon bleed."""

    def test_pentose_not_classified_as_nonose(self):
        """五元环糖 + 长碳链苷元: 糖碳计数不应超过 5-6"""
        from monosaccharide_identifier import _countFragmentCarbons
        from glycan_topology import find_mapped_sugar_units

        # D-Rib + 长碳链苷元 (模拟类胡萝卜素糖苷) [TEST DATA ONLY]
        # 使用简单的 Ribose + 十碳链
        smiles = "OC[C@H]1OC(OCCCCCCCCCC)[C@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pytest.skip("SMILES parse failed")

        units = find_mapped_sugar_units(mol)
        if not units:
            pytest.skip("No sugar units detected")

        ringAtoms = units[0]['ring_atoms']
        carbonCount = _countFragmentCarbons(mol, ringAtoms)

        # 五元环糖的碳计数: 最多 4 (ring) + 1 (C5-OH) = 5
        # 不应超过 9 (Non 的阈值)
        assert carbonCount < 9, (
            f"Pentose carbon count {carbonCount} should not reach nonose threshold (≥9)"
        )

    def test_hexose_carbon_count_reasonable(self):
        """六元环 Glc + 长碳链: 碳计数不应超过 8"""
        from monosaccharide_identifier import _countFragmentCarbons
        from glycan_topology import find_mapped_sugar_units

        # Glc + C20 链 [TEST DATA ONLY]
        smiles = "OC[C@H]1OC(OCCCCCCCCCCCCCCCCCCCC)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pytest.skip("SMILES parse failed")

        units = find_mapped_sugar_units(mol)
        if not units:
            pytest.skip("No sugar units detected")

        ringAtoms = units[0]['ring_atoms']
        carbonCount = _countFragmentCarbons(mol, ringAtoms)

        # Glc 碳计数: 最多 5 (ring C) + 1 (C6) = 6
        # 绝对不应超过 9
        assert carbonCount < 9, (
            f"Hexose carbon count {carbonCount} should not reach nonose threshold"
        )

    def test_sialic_acid_carbon_count_correct(self):
        """唾液酸 (Neu5Ac, 9C) 应正确计数为 ≤9"""
        from monosaccharide_identifier import _countFragmentCarbons

        # Neu5Ac [TEST DATA ONLY]
        smiles = "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@@H](O)[C@@H](O)CO"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pytest.skip("SMILES parse failed")

        ri = mol.GetRingInfo()
        rings = [
            list(r) for r in ri.AtomRings()
            if len(r) in (5, 6) and
            sum(1 for idx in r if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8) >= 1
        ]
        if not rings:
            pytest.skip("No sugar ring found")

        carbonCount = _countFragmentCarbons(mol, rings[0])
        # Neu5Ac 有 9 个碳 (如果全部在 C-C 连续骨架上)
        # 由于 maxHops=4 和 N 截断, 实际计数应 ≤ 9
        assert carbonCount <= 9, (
            f"Sialic acid carbon count {carbonCount} should not exceed 9"
        )


# =====================================================================
# Bug 5: 阈值一致性 (Threshold Consistency)
# =====================================================================

class TestThresholdConsistency:
    """两套分类系统的阈值应保持一致。
    The two classification systems should use consistent thresholds."""

    def test_max_substituent_size_value(self):
        """classify_sugar_parts 的 MAX_SUBSTITUENT_SIZE 应为 10"""
        # 间接验证: 11 原子分支应归入 aglycone
        from glycan_topology import classify_sugar_parts

        smiles = "OC[C@H]1OC(OCCCCCCCCCCC)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pytest.skip("SMILES parse failed")

        classification = classify_sugar_parts(mol)
        aglyconAtoms = classification['aglycone_atoms']

        assert len(aglyconAtoms) > 0, (
            "11-atom branch should be classified as aglycone (threshold=10)"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
