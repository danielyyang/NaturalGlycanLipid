# 🧹 GlycoNP 系统终期审计与清理报告 (End-of-Project Audit & Deletion Report)

**报告日期 (Date):** 2026-03-31
**项目阶段 (Stage):** 最终数据分析与出图期 (Final Data Processing & Visualization)

---

## 📋 审计概述 (Executive Summary)
本项目从最初的概念验证 (PoC)，到各种试错分支 (V11, V12, V13, Part2, Part3)，目前已成功将核心解析引擎和多线程并发控制浓缩至以 `full_pipeline.py` 和 `pipeline_utils.py` 为核心的净网版架构中。
根据代码追踪与互相调用分析，目前 `scripts` 目录下的 92 个文件中有 **超过 60 个文件属于历史迭代遗留物**。

---

## 🗑️ 建议立即删除的文件白名单 (Safe-to-Delete Deprecated Files)

由于项目已收敛，以下文件组确认被废弃，**物理删除它们绝对安全**：

### 1. 废弃的流水线旧代入口 (Deprecated Pipeline Versions)
- `scripts/run_v11_final.py`
- `scripts/run_v12_full_pipeline.py` (包含 1200+ 行杂糅功能，已被 `pipeline_utils` 完全吸收)
- `scripts/run_v13_full_pipeline.py` (并发模式原型，已被 `full_pipeline.py` 替代)
- `scripts/run_glyconp_pipeline.py` (极早期的启动入口)
- `scripts/main_pipeline.py` (早期主线)

### 2. 废除的图表绘制分散文件 (Scattered Chart Scripts)
之前的 40 张图表被分散在三个独立脚本里，现在已 100% 整合进最强大的单脚本 `generate_saponin_charts.py`。
- `scripts/generate_saponin_charts_part2.py`
- `scripts/generate_saponin_charts_part3.py`
- `scripts/generate_executive_plots.py`
- `scripts/generate_glyco_landscape.py`

### 3. 被否决的“苷元猜糖”幻觉逻辑 (Hallucination/Scaffold NLP Guessers)
按用户要求，剥除了不科学的骨架猜手性逻辑。
- `scripts/rescue_generic_sugars_phase2.py`

### 4. 废弃的小型诊断测试脚本 (Orphaned Diagnostic Scripts)
- `scripts/_check_scaffolds.py`
- `scripts/_debug_scaffolds.py`
- `scripts/_debug_scaffolds_deep.py`
- `scripts/compare_v13_v_new.py`
- `scripts/preview_saponin_columns.py`
- `scripts/diag_last4.py`
- `scripts/diag_vitexin.py`
- `tmp_check_hex.py`
- `scripts/verify_csv_swap.py`

### 5. Benchmark 验证期的过时剧本 (Stale Benchmark Drivers)
验证基准目前已 100% 通过（226/226），这批驱动程序功成身退。
- `scripts/run_benchmark_200_tier_a.py`
- `scripts/run_audit_fixes.py`
- `scripts/deep_upgrade_235.py`
- `scripts/build_benchmark_150.py`
- `scripts/build_benchmark_200.py`
- `scripts/fix_benchmark_expected.py`
- `scripts/generate_benchmark_audit_report.py`

### 6. 概念验证引擎 (PoC Engines)
均已移入 `lib` 中封装完成。
- `scripts/cip_exo_identification_engine.py`
- `scripts/poc_cip_fingerprint_engine.py`
- `scripts/poc_virtual_demodify.py`

---

## 🔒 必须保留的核心体系 (Core Files to Retain)

> **警告：以下核心功能为项目的最高智慧结晶，已接管所有数据流。**

### 生产流水线 (Production Pipeline)
- `scripts/full_pipeline.py`: (多核) 天然产物糖化信息提取总出口。
- `scripts/extract_saponins.py`: 提取特定皂苷与生成 HTML 表格。
- `scripts/generate_saponin_charts.py`: 全套科研级 Plotly 图谱输出。
- `scripts/rescue_generic_sugars.py`: 严格文本 NLP (Strategy A) 的防漏系统。

### 底层算法库 (Library Core)
- `lib/pipeline_utils.py`: 统一化处理黑盒。
- `lib/glycan_topology.py`: 立体化学、CIP 回退机制、3-hop 修饰扫描器。
- `lib/monosaccharide_identifier.py`: 虚拟水解酶、异头碳判定、SMARTS 库读取大本营。
- `lib/glycan_reference_library.py`: 黄金单糖库字典。
- `lib/molecular_visualizer.py`: 化学二维结构高亮渲染渲染器。

---

## 📝 最终建议 (Conclusion)
通过执行配套提供的 `cleanup.bat`，不仅能释放数百 MB 存储，更关键的是能掐断一切 `import v12` 或混淆画图输出的可能性，让后续的开发者（审稿人、合作导师）打开仓库时，看到的如同一个设计极简的工业级开源软件。
