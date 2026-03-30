"""
Stereo-Rescue 独立执行脚本 (Standalone Stereo-Rescue Script)
仅运行 Name-Inferred rescue 模块, 不重新执行 rematch
Only runs the Name-Inferred rescue module, skips rematch
"""
import os
import sys
import time
import re
from collections import Counter

import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

from lib.stereochemistry_rescue import rescueSugarSequence, GENERIC_SUGAR_LABELS


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  Stereo-Rescue Standalone (Name-Inferred Only)")
    print("  立体化学确证 — 仅使用名称推断模块")
    print("=" * 70)

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows ({time.time()-t0:.1f}s)")

    # 初始化新列 (Initialize new columns)
    df["Rescued_Sugar_Sequence"] = ""
    df["Rescue_Method"] = ""

    # 找出含泛指标签的行 (Find rows with generic labels)
    genericMask = df["Sugar_Sequence"].astype(str).apply(
        lambda s: any(label in s for label in GENERIC_SUGAR_LABELS)
    )
    genericIndices = df.index[genericMask]
    print(f"  Rows with generic labels: {len(genericIndices):,}")

    rescueCount = 0
    rescueNameCount = 0

    for idx in tqdm(genericIndices, desc="  Stereo-Rescue", ncols=80):
        smiles = str(df.at[idx, "canonical_smiles"])
        currSeq = str(df.at[idx, "Sugar_Sequence"])
        name = str(df.at[idx, "name"]) if "name" in df.columns else None
        iupacName = str(df.at[idx, "iupac_name"]) if "iupac_name" in df.columns else None

        rescuedSeq, method = rescueSugarSequence(
            smiles, currSeq, name, iupacName
        )

        if method and rescuedSeq != currSeq:
            df.at[idx, "Rescued_Sugar_Sequence"] = rescuedSeq
            df.at[idx, "Rescue_Method"] = method
            rescueCount += 1
            if "Name" in method:
                rescueNameCount += 1

    print(f"\n  Rescued: {rescueCount:,} / {len(genericIndices):,}")
    print(f"    Name-Inferred: {rescueNameCount:,}")

    # 统计 rescue 后的糖名分布 (Stats of rescued sugar names)
    rescuedSeqs = df["Rescued_Sugar_Sequence"][df["Rescued_Sugar_Sequence"] != ""]
    if len(rescuedSeqs) > 0:
        rescuedTokens = []
        for seq in rescuedSeqs:
            rescuedTokens.extend(re.findall(
                r'Neu5Ac|Neu5Gc|KDO|Hept|Oct|Non|'
                r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*|Hex|dHex|Pen|HexA|HexN|'
                r'HexNAc',
                str(seq)))
        rescuedCounter = Counter(rescuedTokens)
        print(f"\n  Rescued token distribution:")
        for token, count in rescuedCounter.most_common(15):
            print(f"    {token:<20s} {count:>6,}")

    # 保存 (Save)
    outPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")
    print(f"\n  Updated: {outPath}")
    print(f"  Time: {time.time()-t0:.0f}s")
    print(f"  {'='*60}")


if __name__ == "__main__":
    main()
