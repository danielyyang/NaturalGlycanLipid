"""
测试本地 LOTUS 分类学匹配引擎 — 使用 Mock 数据
Test Local LOTUS Taxonomy Matching Engine — Using Mock Data

Mock 数据格式说明 (Mock Data Format):
LOTUS dump CSV 需要以下列 (LOTUS dump CSV requires these columns):
  - structure_inchikey:          完整 InChIKey (27 字符, e.g. "GZUZMTKPFCAKIR-HPODLMMYSA-N")
  - organism_name:               物种名 (e.g. "Panax ginseng")
  - organism_taxonomy_01domain:  域 (e.g. "Eukaryota")
  - organism_taxonomy_02kingdom: 界 (e.g. "Plantae")
  - organism_taxonomy_05family:  科 (e.g. "Araliaceae")
  - organism_taxonomy_06genus:   属 (e.g. "Panax")
  - organism_taxonomy_08species: 种 (e.g. "Panax ginseng")

目标 CSV 需要以下列 (Target CSV requires these columns):
  - standard_inchi_key:          COCONUT 中的 InChIKey
  - organisms:                   现有 organism 值 (可为空)
  - Family:                      现有 Family 值 (可为空)
"""
import os
import sys
import tempfile
import pytest
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.taxonomy_lotus_matcher import loadLotusDump, fillTaxonomyFromLotus


# =====================================================================
# [TEST DATA ONLY] — 合成测试数据
# =====================================================================

def createMockLotusCsv(tmpDir: str) -> str:
    """创建模拟 LOTUS dump CSV / Create mock LOTUS dump CSV"""
    data = {
        "structure_inchikey": [
            # block-1: GZUZMTKPFCAK -> Panax ginseng
            "GZUZMTKPFCAKIR-HPODLMMYSA-N",
            # 同一个 Block-1 的不同立体异构体 (Same Block-1, different stereoisomer)
            "GZUZMTKPFCAKIR-XXXXXXXYSA-N",
            # block-1: REFJWTPEDVJJIY -> quercetin
            "REFJWTPEDVJJIY-UHFFFAOYSA-N",
            # block-1: OVSQVDMCBVZWGM -> kaempferol
            "OVSQVDMCBVZWGM-UHFFFAOYSA-N",
            # block-1: ABCDEFGHIJKLMN -> synthetic compound (no match expected)
            "TESTINCHIKEY01-UHFFFAOYSA-N",
        ],
        "organism_name": [
            "Panax ginseng",
            "Panax notoginseng",
            "Allium cepa",
            "Ginkgo biloba",
            "Streptomyces griseus",
        ],
        "organism_taxonomy_01domain": [
            "Eukaryota", "Eukaryota", "Eukaryota", "Eukaryota", "Bacteria",
        ],
        "organism_taxonomy_02kingdom": [
            "Plantae", "Plantae", "Plantae", "Plantae", "Bacteria",
        ],
        "organism_taxonomy_05family": [
            "Araliaceae", "Araliaceae", "Amaryllidaceae", "Ginkgoaceae", "Streptomycetaceae",
        ],
        "organism_taxonomy_06genus": [
            "Panax", "Panax", "Allium", "Ginkgo", "Streptomyces",
        ],
        "organism_taxonomy_08species": [
            "Panax ginseng", "Panax notoginseng", "Allium cepa", "Ginkgo biloba", "Streptomyces griseus",
        ],
    }
    path = os.path.join(tmpDir, "mock_lotus.csv")
    pd.DataFrame(data).to_csv(path, index=False)
    return path


def createMockTargetDf() -> pd.DataFrame:
    """创建模拟目标 DataFrame (COCONUT 格式) / Create mock target DataFrame"""
    return pd.DataFrame({
        "standard_inchi_key": [
            # 应匹配 Panax ginseng (Match expected: Panax ginseng/notoginseng)
            "GZUZMTKPFCAKIR-DIFFERENTSA-N",
            # 应匹配 quercetin → Allium cepa (Match expected)
            "REFJWTPEDVJJIY-UHFFFAOYSA-N",
            # 无匹配 (No match: InChIKey not in LOTUS)
            "NOMATCHFOUND00-UHFFFAOYSA-N",
            # 已有 organism，不应被覆盖 (Already has organism, should NOT be overwritten)
            "OVSQVDMCBVZWGM-UHFFFAOYSA-N",
        ],
        "organisms": [
            "",         # 空 → 应被填补
            "NULL",     # NULL → 应被填补
            "",         # 空 → 无匹配，保持空
            "Existing Org",  # 已有 → 不覆盖
        ],
        "Family": [
            "",
            "",
            "",
            "ExistingFamily",
        ],
        "name": ["Ginsenoside Rb1", "Quercetin", "Unknown Compound", "Kaempferol"],
    })


# =====================================================================
# Tests
# =====================================================================

class TestLoadLotusDump:
    """测试 LOTUS dump 加载和索引构建 / Test LOTUS dump loading & index building"""

    def test_basic_load(self, tmp_path):
        """基本加载：正确读取并构建 Block-1 索引"""
        lotusCsv = createMockLotusCsv(str(tmp_path))
        index = loadLotusDump(lotusCsv)

        assert isinstance(index, pd.DataFrame)
        assert "inchikey_block1" in index.columns
        assert "organism_name" in index.columns
        # 5 条记录，但 Block-1 GZUZMTKPFCAK 有 2 条 → 合并后应有 4 个唯一 Block-1
        assert len(index) == 4, f"Expected 4 unique Block-1 entries, got {len(index)}"

    def test_aggregation(self, tmp_path):
        """聚合验证：同一 Block-1 的多个 organism 应被 | 连接"""
        lotusCsv = createMockLotusCsv(str(tmp_path))
        index = loadLotusDump(lotusCsv)

        # Block-1 = GZUZMTKPFCAK 应该有 Panax ginseng|Panax notoginseng
        ginsengRow = index[index["inchikey_block1"] == "GZUZMTKPFCAKIR"]
        assert len(ginsengRow) == 1
        orgValue = ginsengRow.iloc[0]["organism_name"]
        assert "Panax ginseng" in orgValue
        assert "Panax notoginseng" in orgValue

    def test_block1_length(self, tmp_path):
        """所有 Block-1 键应为 14 字符"""
        lotusCsv = createMockLotusCsv(str(tmp_path))
        index = loadLotusDump(lotusCsv)

        for val in index["inchikey_block1"]:
            assert len(val) == 14, f"Block-1 should be 14 chars, got {len(val)}: {val}"


class TestFillTaxonomy:
    """测试核心匹配逻辑 / Test core matching logic"""

    @pytest.fixture
    def lotusIndex(self, tmp_path):
        lotusCsv = createMockLotusCsv(str(tmp_path))
        return loadLotusDump(lotusCsv)

    def test_fill_missing_organism(self, lotusIndex):
        """空 organism 应被填补 / Missing organism should be filled"""
        targetDf = createMockTargetDf()
        resultDf, imputedCells = fillTaxonomyFromLotus(targetDf, lotusIndex)

        # Row 0: 空 → 应被填补为 Panax ginseng|Panax notoginseng
        row0Org = str(resultDf.iloc[0]["organisms"])
        assert "Panax" in row0Org, f"Row 0 should have Panax, got: {row0Org}"

    def test_fill_null_organism(self, lotusIndex):
        """'NULL' 标记的 organism 应被填补"""
        targetDf = createMockTargetDf()
        resultDf, _ = fillTaxonomyFromLotus(targetDf, lotusIndex)

        # Row 1: NULL → 应被填补为 Allium cepa
        row1Org = str(resultDf.iloc[1]["organisms"])
        assert "Allium" in row1Org, f"Row 1 should have Allium, got: {row1Org}"

    def test_no_match_stays_empty(self, lotusIndex):
        """无匹配的行保持原值 / No-match rows should stay unchanged"""
        targetDf = createMockTargetDf()
        resultDf, _ = fillTaxonomyFromLotus(targetDf, lotusIndex)

        # Row 2: 无匹配 → 保持空
        row2Org = str(resultDf.iloc[2]["organisms"])
        assert row2Org in ("", "nan", "NULL"), f"Row 2 should be empty, got: {row2Org}"

    def test_existing_organism_not_overwritten(self, lotusIndex):
        """已有 organism 的行不应被覆盖 / Existing organism should NOT be overwritten"""
        targetDf = createMockTargetDf()
        resultDf, _ = fillTaxonomyFromLotus(targetDf, lotusIndex)

        # Row 3: 已有 "Existing Org" → 不覆盖
        row3Org = str(resultDf.iloc[3]["organisms"])
        assert row3Org == "Existing Org", f"Row 3 should keep 'Existing Org', got: {row3Org}"

    def test_family_filled(self, lotusIndex):
        """空 Family 应被填补 / Missing Family should be filled"""
        targetDf = createMockTargetDf()
        resultDf, _ = fillTaxonomyFromLotus(targetDf, lotusIndex)

        # Row 0: 应被填补为 Araliaceae
        row0Fam = str(resultDf.iloc[0]["Family"])
        assert "Araliaceae" in row0Fam, f"Row 0 Family should have Araliaceae, got: {row0Fam}"

    def test_existing_family_not_overwritten(self, lotusIndex):
        """已有 Family 不应被覆盖 / Existing Family should NOT be overwritten"""
        targetDf = createMockTargetDf()
        resultDf, _ = fillTaxonomyFromLotus(targetDf, lotusIndex)

        row3Fam = str(resultDf.iloc[3]["Family"])
        assert row3Fam == "ExistingFamily", f"Row 3 Family should keep 'ExistingFamily', got: {row3Fam}"

    def test_imputed_cells_tracking(self, lotusIndex):
        """被填补的单元格应被正确追踪 / Imputed cells should be tracked"""
        targetDf = createMockTargetDf()
        _, imputedCells = fillTaxonomyFromLotus(targetDf, lotusIndex)

        # 至少 Row 0 和 Row 1 的 organisms 应被标记
        orgImputed = [c for c in imputedCells if c[1] == "organisms"]
        assert len(orgImputed) >= 2, f"Expected >= 2 organism imputations, got {len(orgImputed)}"

    def test_no_for_loop_performance(self, lotusIndex):
        """性能验证：大数据量下应在 1 秒内完成
        Performance: should complete in < 1 second for large datasets"""
        import time
        # 创建 10,000 行测试数据 (Create 10K row test data)
        largeDf = pd.DataFrame({
            "standard_inchi_key": [f"GZUZMTKPFCAKIR-{i:010d}SA-N" for i in range(10000)],
            "organisms": [""] * 10000,
            "Family": [""] * 10000,
        })

        startTime = time.time()
        resultDf, _ = fillTaxonomyFromLotus(largeDf, lotusIndex)
        elapsed = time.time() - startTime

        assert elapsed < 2.0, f"10K rows should complete in < 2s, took {elapsed:.2f}s"
        assert len(resultDf) == 10000


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
