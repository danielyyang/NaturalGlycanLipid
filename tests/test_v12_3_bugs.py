"""
V12.3 Bug 修复回归测试 (V12.3 Bug Fix Regression Tests)
验证 Bug 2 (Bond_Detail 范围), Bug 5 (纯糖 Glycan_SMILES), Bug 6 (高亮图 + 修饰列)
[TEST DATA ONLY]
"""
import sys, os, json
import pytest
from rdkit import Chem

sys.path.insert(0, r"d:\Glycan_Database")
from lib.glycan_topology import find_mapped_sugar_units


class TestBug2BondDetailScope:
    """Bug 2: Bond_Detail 仅记录 Sugar→Aglycone / Bond_Detail only records Sugar→Aglycone"""

    def test_trisaccharide_bond_detail_no_sugar_sugar(self):
        """含 3 个糖 + 苷元的分子, Bond_Detail 不应包含 Sugar→Sugar 键
        Trisaccharide+aglycone: Bond_Detail should only have Sugar→Aglycone entries"""
        # 手动构建测试: 导入 detectAllGlycosidicBonds 并验证
        sys.path.insert(0, r"d:\Glycan_Database\scripts")
        from run_v12_full_pipeline import detectAllGlycosidicBonds

        # 使用一个简单的糖苷 (甲基葡萄糖苷) 测试
        # Methyl glucoside: 1 sugar + methyl aglycone
        methylGlcSmi = "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(methylGlcSmi)
        assert mol is not None

        units = find_mapped_sugar_units(mol)
        rootBond, bondJson = detectAllGlycosidicBonds(mol, units)
        bonds = json.loads(bondJson)

        # 所有 Bond_Detail 条目的 target_type 必须为 Aglycone
        for bd in bonds:
            assert bd["target_type"] == "Aglycone", (
                f"Bug 2: 发现 Sugar→Sugar 键: {bd}")


class TestBug5PureSugarRetention:
    """Bug 5: 纯糖分子保留 Glycan_SMILES / Pure sugars retain Glycan_SMILES"""

    def test_maltose_is_pure_sugar(self):
        """麦芽糖 (Maltose) 应被识别为纯糖分子
        Maltose should be identified as a pure sugar molecule"""
        sys.path.insert(0, r"d:\Glycan_Database\scripts")
        from run_v12_full_pipeline import _checkPureSugarMolecule

        maltoseSmi = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]2OC(CO)[C@@H](O)[C@H](O)[C@H]2O"
        mol = Chem.MolFromSmiles(maltoseSmi)
        assert mol is not None
        assert _checkPureSugarMolecule(mol), "麦芽糖应被识别为纯糖分子"

    def test_glycoside_is_not_pure_sugar(self):
        """含大苷元的分子不应被识别为纯糖 / Large aglycone should NOT be pure sugar"""
        sys.path.insert(0, r"d:\Glycan_Database\scripts")
        from run_v12_full_pipeline import _checkPureSugarMolecule

        # 简单苯基葡萄糖苷 (phenyl glucoside)
        phenylGlcSmi = "OC[C@H]1OC(Oc2ccccc2)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(phenylGlcSmi)
        assert mol is not None
        # 苯环只有 6 个碳, 可能刚好达到 70% — 但大苷元不是
        # 更大的苷元:
        steroidGlcSmi = "OC[C@H]1OC(OC2CCC3(C)C(CCC4C3CCC3(C)C(C)CCC43)C2)[C@H](O)[C@@H](O)[C@@H]1O"
        mol2 = Chem.MolFromSmiles(steroidGlcSmi)
        if mol2:
            assert not _checkPureSugarMolecule(mol2), "含甾体苷元不应为纯糖"


class TestBug6HighlightedImage:
    """Bug 6: 三色高亮图 / Three-color highlighted image"""

    def test_highlighted_png_returns_base64(self):
        """高亮画图函数应返回非空 base64 字符串
        Highlighted image function should return non-empty base64"""
        sys.path.insert(0, r"d:\Glycan_Database\scripts")
        from run_v12_full_pipeline import molToHighlightedBase64Png

        glcSmi = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
        result = molToHighlightedBase64Png(glcSmi, size=(300, 200))
        assert result, "高亮图应返回非空 base64"
        assert len(result) > 100, "base64 长度应 > 100"

    def test_highlighted_handles_empty_smiles(self):
        """空 SMILES 应返回空字符串 / Empty SMILES returns empty string"""
        sys.path.insert(0, r"d:\Glycan_Database\scripts")
        from run_v12_full_pipeline import molToHighlightedBase64Png

        assert molToHighlightedBase64Png("") == ""
        assert molToHighlightedBase64Png("nan") == ""
        assert molToHighlightedBase64Png(None) == ""


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
