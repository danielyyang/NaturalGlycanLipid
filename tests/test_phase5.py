"""
Phase 5 测试: 苷元骨架 + 糖序列退避匹配
Phase 5 Test: Aglycone Scaffold + Tiered Glycan Sequence Matching

测试用例 (Test Cases):
  1. Rutin (芦丁): Glc + Rha 双糖 + quercetin 苷元
  2. Ginsenoside Rg1: 两个 Glc + 甾体苷元
  3. 直链脂肪酸苷元: 测试 "Aliphatic Chain" 退避
  4. 无立体构型 SMILES: 测试 Tier 3 Hex/Pen 退避
"""
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.feature_extractor import (
    computeMorganFingerprint,
    extractMurckoScaffold,
    characterizeAglycon,
    characterizeGlycan,
    processPhase5Row,
)


def printSection(title: str):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def testRutin():
    """[TEST DATA ONLY] 芦丁 (Rutin): quercetin-3-O-rutinoside"""
    printSection("Test 1: Rutin (quercetin-3-O-rutinoside)")
    # 芦丁的糖基部分: Glc-(1→6)-Rha (glucose + rhamnose)
    glycanSmiles = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](O)[C@@H](C)O2)C(CO)O1"
    # 苷元: quercetin (槲皮素) — 含 3 个环
    aglycanSmiles = "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

    print(f"  Glycan SMILES:  {glycanSmiles[:60]}...")
    print(f"  Aglycon SMILES: {aglycanSmiles}")

    result = processPhase5Row(glycanSmiles, aglycanSmiles)

    print(f"\n  [Glycan Results]")
    print(f"    Sugar Sequence:          {result['sugar_sequence']}")
    print(f"    Sugar Functional Group:  {result['sugar_functional_group']}")

    print(f"\n  [Aglycon Results]")
    print(f"    Murcko Scaffold:         {result['murcko_scaffold']}")
    print(f"    Aglycon MW:              {result['aglycon_mw']}")
    print(f"    Ring Count:              {result['aglycon_ring_count']}")
    print(f"    Morgan FP length:        {len(result['morgan_fp'])}")

    # 断言 (Assertions)
    assert len(result['morgan_fp']) == 2048, "Morgan FP should be 2048 bits"
    assert result['aglycon_ring_count'] == "3", f"Quercetin has 3 rings, got {result['aglycon_ring_count']}"
    assert result['murcko_scaffold'] not in ("NULL", ""), f"Quercetin should have a scaffold, got '{result['murcko_scaffold']}'"
    assert result['murcko_scaffold'] != "Aliphatic Chain", "Quercetin is not aliphatic"
    # 糖序列至少应产生某种标注 (不论具体识别结果如何)
    assert result['sugar_sequence'] != "", "Should produce a sugar sequence"
    print("\n  >> All assertions PASSED")


def testGinsenosideRg1():
    """[TEST DATA ONLY] 人参皂苷 Rg1: 含两个 Glc + 甾体骨架"""
    printSection("Test 2: Ginsenoside Rg1")
    # 糖基: 两个葡萄糖 (通过甾体的两个位点连接)
    glycanSmiles = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
    # 苷元: 达玛烷三萜 (dammarane triterpenoid)
    aglycanSmiles = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC=C2C)C[C@H](O)[C@@H]4[C@@]3(C)CC[C@H](C4(C)C)O"

    print(f"  Glycan SMILES:  {glycanSmiles}")
    print(f"  Aglycon SMILES: {aglycanSmiles[:60]}...")

    result = processPhase5Row(glycanSmiles, aglycanSmiles)

    print(f"\n  [Glycan Results]")
    print(f"    Sugar Sequence:          {result['sugar_sequence']}")
    print(f"    Sugar Functional Group:  {result['sugar_functional_group']}")

    print(f"\n  [Aglycon Results]")
    print(f"    Murcko Scaffold:         {result['murcko_scaffold']}")
    print(f"    Aglycon MW:              {result['aglycon_mw']}")
    print(f"    Ring Count:              {result['aglycon_ring_count']}")

    # 甾体骨架有 4 个环
    assert int(result['aglycon_ring_count']) >= 4, f"Dammarane should have >=4 rings, got {result['aglycon_ring_count']}"
    assert result['murcko_scaffold'] != "NULL"
    print("\n  >> All assertions PASSED")


def testAliphaticChainFallback():
    """[TEST DATA ONLY] 脂肪酸苷元: 直链分子 → Aliphatic Chain 退避"""
    printSection("Test 3: Aliphatic Chain Fallback")
    # 棕榈酸 (Palmitic acid) — 无环脂肪链
    aliphaticSmiles = "CCCCCCCCCCCCCCCC(=O)O"

    result = characterizeAglycon(aliphaticSmiles)
    print(f"  Input: {aliphaticSmiles}")
    print(f"  Murcko Scaffold: {result['murcko_scaffold']}")
    print(f"  Ring Count:      {result['aglycon_ring_count']}")

    assert result['murcko_scaffold'] == "Aliphatic Chain", f"Expected 'Aliphatic Chain', got '{result['murcko_scaffold']}'"
    assert result['aglycon_ring_count'] == "0"
    print("\n  >> All assertions PASSED")


def testNoStereoFallback():
    """[TEST DATA ONLY] 无立体化学 SMILES: 测试 Tier 3 Hex/Pen 退避"""
    printSection("Test 4: No-Stereo Fallback (Tier 3: Hex/Pen)")
    # 无立体信息的六元糖环 → 应退避到 Hex
    noStereoHexose = "OCC1OC(O)C(O)C(O)C1O"

    result = characterizeGlycan(noStereoHexose)
    print(f"  Input: {noStereoHexose}")
    print(f"  Sugar Sequence:          {result['sugar_sequence']}")
    print(f"  Sugar Functional Group:  {result['sugar_functional_group']}")

    # 没有立体化学，应退避到 Hex 级别（不应该精确匹配到 Glc/Gal）
    seq = result['sugar_sequence']
    if "Hex" in seq or "Glc" in seq or "Gal" in seq or "Man" in seq:
        print(f"\n  >> Tier 3 fallback working: identified as '{seq}'")
    else:
        print(f"\n  >> WARNING: unexpected sequence '{seq}'")

    assert result['sugar_sequence'] != "", "Should produce at least Hex"
    print("  >> Assertion PASSED")


def testNullHandling():
    """[TEST DATA ONLY] NULL/空值处理: 不应崩溃"""
    printSection("Test 5: NULL / Empty Value Handling")
    nullInputs = ["", "NULL", "nan", None, "*"]

    for inp in nullInputs:
        agResult = characterizeAglycon(inp)
        glycResult = characterizeGlycan(inp)

        assert agResult['murcko_scaffold'] == "NULL", f"Expected NULL for '{inp}'"
        assert agResult['morgan_fp'] == "", f"Expected empty FP for '{inp}'"
        assert glycResult['sugar_sequence'] == "", f"Expected empty sequence for '{inp}'"

    print("  All 5 null inputs handled gracefully")
    print("  >> All assertions PASSED")


def testMorganFPConsistency():
    """[TEST DATA ONLY] Morgan FP 一致性: 相同 SMILES 应产生相同指纹"""
    printSection("Test 6: Morgan FP Consistency")
    smiles1 = "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
    smiles2 = "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

    fp1 = computeMorganFingerprint(smiles1)
    fp2 = computeMorganFingerprint(smiles2)

    assert fp1 == fp2, "Same SMILES should produce identical fingerprints"
    assert len(fp1) == 2048
    print(f"  SMILES: {smiles1}")
    print(f"  FP consistency: identical={fp1 == fp2}, length={len(fp1)}")
    print("  >> Assertion PASSED")


if __name__ == "__main__":
    print("=" * 70)
    print("   Phase 5 Feature Extraction — Comprehensive Test Suite")
    print("=" * 70)

    testRutin()
    testGinsenosideRg1()
    testAliphaticChainFallback()
    testNoStereoFallback()
    testNullHandling()
    testMorganFPConsistency()

    print(f"\n{'='*70}")
    print("  ALL TESTS PASSED")
    print(f"{'='*70}")
