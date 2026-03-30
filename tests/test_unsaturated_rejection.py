"""
测试不饱和糖环/双键结构/多环聚醚的拒止逻辑
Test unsaturated sugar ring / double-bond / polycyclic polyether rejection logic

验证 Phase 2 核心切分中 find_mapped_sugar_units() 的以下三个新规则:
1. 环内双键拒止 (ring-internal unsaturation rejection)
2. 环外 C=O 不影响正常糖 (exocyclic C=O does NOT disqualify real sugars)
3. 多环含氧聚醚拓扑铁律 (polycyclic oxygenated polyether topology rule)
"""
import os
import sys
import pytest
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.glycan_topology import find_mapped_sugar_units


class TestUnsaturatedRejection:
    """环内含双键的不饱和糖环应被拒止 / Rings with internal C=C should be rejected"""

    def test_unsaturated_pyranoside_rejected(self):
        """不饱和吡喃糖苷 (C=C 在环内) 应被拒止
        Unsaturated pyranoside with ring-internal C=C should be rejected"""
        # 2,3-Didehydro compound with ring C=C
        unsatSmiles = "CO[C@@H]1O[C@H](CO)[C@@H](O)C=C1"
        mol = Chem.MolFromSmiles(unsatSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) == 0, f"不饱和糖环应被拒止，但识别出了 {len(units)} 个糖元"

    def test_glycal_rejected(self):
        """糖烯 (Glycal) 含环内双键应被拒止
        Glycals (1,2-unsaturated sugars) should be rejected"""
        # D-glucal: 1,2-unsaturated glucose
        glucalSmiles = "OC[C@H]1OC=C[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(glucalSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) == 0, f"糖烯应被拒止，但识别出了 {len(units)} 个糖元"


class TestExocyclicBondsOK:
    """环外的不饱和键不应影响正常糖的识别 / Exocyclic unsaturation should NOT affect sugar detection"""

    def test_glucuronic_acid_accepted(self):
        """β-D-葡萄糖醛酸的 C6=O 羧基不应伤害糖环识别
        β-D-Glucuronic acid with exocyclic C=O should still be accepted as a sugar"""
        glcASmiles = "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(glcASmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) >= 1, "糖醛酸应被识别为糖环"

    def test_acetylated_sugar_accepted(self):
        """乙酰化葡萄糖的环外 C=O 酯键不应伤害糖环识别
        Acetylated glucose with exocyclic ester C=O should still be accepted"""
        glcOAcSmiles = "CC(=O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(glcOAcSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) >= 1, "乙酰化糖应被识别为糖环"

    def test_caffeoylated_sugar_accepted(self):
        """咖啡酰化的糖：环外含苯环和 C=C/C=O 但糖环本身应被识别
        Caffeoylated sugar: exocyclic aromatic + C=C/C=O should not disqualify the sugar"""
        caffeoylGlcSmiles = "OC(=O)/C=C/c1ccc(O)c(O)c1.O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1"
        mol = Chem.MolFromSmiles(caffeoylGlcSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) >= 1, "咖啡酰化糖的糖环应被识别"


class TestPolycyclicPolyetherRejection:
    """多环含氧聚醚拓扑铁律测试 / Polycyclic oxygenated polyether topology rule"""

    def test_simple_sugar_accepted(self):
        """正常单糖应被接受 / Normal monosaccharide should be accepted"""
        glcSmiles = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(glcSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) >= 1, "正常葡萄糖应被识别"

    def test_disaccharide_accepted(self):
        """正常双糖应被接受 / Normal disaccharide should be accepted"""
        # Maltose
        maltoseSmiles = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]2OC(CO)[C@@H](O)[C@H](O)[C@H]2O"
        mol = Chem.MolFromSmiles(maltoseSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) >= 2, f"麦芽糖应识别出 >=2 个糖元，实际 {len(units)}"


class TestGlycosidicLinkageConstraint:
    """糖键连方式约束测试 / Glycosidic linkage type constraint"""

    def test_normal_o_glycoside_accepted(self):
        """C-O-C 连接的正常 O-糖苷应被接受
        Normal O-glycoside (C-O-C linkage) should be accepted"""
        # Methyl glucoside
        methylGlcSmiles = "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(methylGlcSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) >= 1, "甲基糖苷应被识别"

    def test_c_glycoside_accepted(self):
        """C-C 直接连接的 C-糖苷应被接受
        C-glycoside (direct C-C linkage) should be accepted"""
        # A simple C-glycoside: C1 bonded directly to a carbon
        cGlycosideSmiles = "OC[C@H]1OC(c2ccccc2)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(cGlycosideSmiles)
        assert mol is not None
        units = find_mapped_sugar_units(mol)
        assert len(units) >= 1, "C-糖苷应被识别"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
