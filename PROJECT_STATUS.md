# 📈 GlycoNP 项目阶段状态报告 (Project Status)

**更新日期 (Date Updated): 2026-03-31**
**当前阶元 (Current Phase): 终点数据汇编与可视化输出阶元 (Final Data Assembly & Analytics Publishing)**

---

## 🏃 进展速览 (At a Glance)
本项目历次大重构及试错版本探索（V11, V12, V13）已经圆满收结。目前主干执行框架（以 `full_pipeline.py` 与 `pipeline_utils.py` 为骨架的并行分析引擎）已进入封卷阶段。
- **全量化学拆分与特征比对已完工。**
- **Saponin 子集的定制化分析表和高维统计图谱已经全部自动化。**
- **代码库已经过全面的大扫除与历史版本清理（Audit Report 2026-03-31）。**

### 🏆 研发成果 (Key Deliverables Reached)
1. **多核超大并发处理 (Concurrency Enabled):** V13 理念被合并进最终版，管线耗时从单进程 2 小时压制成近线性的时间消耗（~5 分钟）。
2. **纯粹的 CIP/拓扑立体分析 (Stereocenter Integrity):** 去除了依赖分类学骨架来盲推糖环手性的幻觉（Hallucination）风险。目前的 ~15 万 `D-Glc` 基团，有近乎 >99.7% 的绝对置信度来源于 RDKit `[C@@H]` 基元识别与明确的文献名称文本提取。
3. **彻底打通 3-hop Scanner 硫酸化识别**: 能够精确捕集 C6 甚至脱氧糖 C4 尾端附着的外延酯/磷酸、氨基，完美地重构了脱元后的化学全骨架。

---

## 📅 下一步行动 (Next Steps: Paper Preparation)

由于计算流程已定格，下面的进展应完全围绕文章/报告的撰写与出版：

- [ ] **Data Publishing 验证**:
  梳理提取的 `GlycoNP_Saponin_DB.csv` 里的新颖皂苷序列，与学术界过往认知的皂苷连接规律进行对比描述，发掘反常理或极具药学价值的新颖糖修饰 (Rare Modification & Branched Co-occurrence)。
- [ ] **Figure Selection**:
  管线自动化产出的 40 张科学图谱中，挑选核心趋势。比如：A1 的极纯净单糖排名（摒弃了 3300 多个骨架幻觉），S10 的共发生网络图，呈现天然产物对修饰物种选择的高度偏好。
- [ ] **Methods Section 编写**: 
  参考 `lib` 下的方法名、`SizeThreshold=3` 的异头隔离容忍度设计，以及 `rescue_generic_sugars.py` 中唯二依仗的文本分析 NLP 兜底规则，作为算法方法学写入补充材料 (SI)。

## 🚨 已处理的技术负债 (Resolved Technical Debt)

* **版本混乱 (Version Sprawl):** 所有图表、输入输出脚本全部断开了对诸如 `_v12`, `_V13_Pruned` 的依赖硬编码。现在接口固定为 `GlycoNP_Deep_Enriched_Final.csv`。
* **伪精确数据注入 (False-Positive NLP Guesses):** V12-V13 历史阶段用 `Strategy C` (骨架特征统计关联猜测) 的代码已经被**彻 底 切 除**。这捍卫了项目作为数据标定平台的底层科学伦理，保留不足 1% 比例的 2D `Hex` 本为真理面貌的直接呈现。
* **文件冗余积压 (File Hoarding):** 几十个探索用 `_diag_` 的 Python 测试用例与作废流水线已被清出当前 GitHub 索引目录，释放存储的同时提供了即时读图的便利性。
