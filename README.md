# GlycoNP-Pipeline

## 这个项目在做什么？ (What Does This Project Do?)

**GlycoNP-Pipeline** 是一个从天然产物数据库 (COCONUT) 中系统性挖掘**糖缀合物** (Glycoconjugates) 规律的计算化学管线。

在天然产物化学中，糖链不是装饰品 —— 它决定了化合物的水溶性、细胞识别、生物利用度甚至药理活性。但传统的天然产物研究几乎都聚焦在苷元骨架上，**糖链被当作噪声丢弃了**。

我们的目标是：把这些被忽视的糖链"翻出来"，看看大自然在不同物种、不同骨架、不同生物功能的化合物上，**到底选择了哪些糖、做了哪些修饰、为什么**。

## What question are we answering?

> **"在 94,242 个天然含糖化合物中，特定的糖链序列和修饰模式，是否与物种来源、苷元骨架类型存在系统性的关联？"**

---

## 七阶段处理流程 (Seven-Phase Pipeline)

```
原始数据 ──► Phase 1 ──► Phase 2 ──► Phase 3 ──► Phase 4
                                                     │
  报告  ◄── Phase 7 ◄── Phase 6 ◄── Phase 5 ◄───────┘
```

### Phase 1: 含糖筛选 (Sugar Filtering) ✅
从 COCONUT 全量数据库中筛选含糖化合物。基于 RDKit 糖环检测 (五元/六元含氧杂环 + 多羟基特征)，筛出 94,242 条，按 `standard_inchi_key` 严格去重。

### Phase 2: 糖-苷元精确切分 (Glycan–Aglycon Cleavage) ✅
异头碳定位 → 糖苷键识别 → BFS 可达性判定 → 碳原子守恒断言。纯糖分子检测 (≥70% 重原子属于糖域时不切分)。

### Phase 3: 核苷酸糖 & 糖肽识别 (Nucleotide Sugar & Glycopeptide Detection) ✅
SMARTS 匹配碱基 (嘌呤/嘧啶) + 磷酸基团共存检测。肽键明确排除 NAc (N-乙酰基) 甲基碳。

### Phase 4: 分类学填补 (Taxonomy Enrichment) ⚠️
LOTUS InChIKey 映射 + PubChem/GBIF 在线查询 + 240+ 植物属/50+ 真菌属/36+ 细菌属字典。当前 ~49.4% 化合物缺少 organism 信息。

### Phase 5: 糖序列 & 苷元骨架提取 (Sugar Sequence & Scaffold Extraction) ✅
**三层匹配引擎** (Tier 1: 绝对手性 → Tier 2: 异头碳特赦 → Tier 3: 模糊手性) + **CIP+Exo 指纹鉴定引擎** v2.0 + 虚拟脱修饰 + NLP 名称救援。120+ 参考糖全库。

### Phase 6: 智能化学分类 (Intelligent Chemical Classification) ✅
SMARTS 骨架规则 + 糖脂特异捕获 + Tanimoto 最近邻回退。

### Phase 7: 可视化与报告 (Visualization & Reports) ✅
三色高亮体系 (🔴糖环 / 🟡修饰 / 🔵苷元) + HTML/Excel 报告 + 糖苷键 α/β 标注。

---

## 核心引擎架构 (Core Engine Architecture)

```
                    ┌──────────────────────────────────────┐
                    │   lib/monosaccharide_identifier.py    │
                    │   ── 单糖鉴定主引擎 v10 ──           │
                    │   generate_refined_sequence()         │
                    │    └→ identify_monosaccharide_v10()   │
                    │        └→ identify_monosaccharide_v2()│
                    └───────────────┬──────────────────────┘
                                    │
                ┌───────────────────┴──────────────────────┐
                │                                          │
    ┌───────────▼───────────┐          ┌───────────────────▼──────┐
    │ lib/cip_exo_engine.py │          │ lib/virtual_demodify.py  │
    │ CIP+Exo 指纹鉴定      │─────────→│ 虚拟脱修饰 (19 SMIRKS)   │
    │ (walkSugarRing +      │          │ + 3 保护门               │
    │  extractFingerprint)  │          └──────────────────────────┘
    └───────────────────────┘
                │
    ┌───────────▼────────────────┐
    │ lib/glycan_reference_      │
    │    library.py              │
    │ 120+ 参考糖 SMILES (KEGG)  │
    └────────────────────────────┘
```

### 公开 API (调用优先级)

| API | 用途 | 调用场景 |
|-----|------|----------|
| `generate_refined_sequence(mol)` | 完整糖序列 + 修饰格式化 | **生产管线** (V12/V13) |
| `analyze_glycan(smiles)` | SMILES → 序列 (便捷包装) | 基准测试 / Demo |
| `identify_monosaccharide_v10(mol, ring)` | 单环鉴定 (含 CIP/Exo + Rescue) | 调试 / 基准测试逐环分析 |
| `identify_monosaccharide_v2(mol, ring)` | 基础匹配 (仅 SMARTS) | **仅供 v10 内部调用** |

> ⚠️ **约定**: 所有新脚本必须使用 `generate_refined_sequence()` 或 `analyze_glycan()`。禁止直接调用 `v2`。

---

## 项目目录结构 (Project Structure)

```
D:\Glycan_Database\
├── lib/                               # 🟢 核心算法层 (15 模块)
│   ├── glycan_topology.py             #   糖环检测 + 拓扑验证 + 切分 (58KB)
│   ├── monosaccharide_identifier.py   #   三层匹配 + CIP/Exo + Rescue (106KB)
│   ├── glycan_reference_library.py    #   120+ 参考糖 SMILES/SMARTS (37KB)
│   ├── cip_exo_engine.py             #   CIP+Exo 指纹鉴定引擎 (27KB)
│   ├── virtual_demodify.py           #   虚拟脱修饰 + 保护门 (17KB)
│   ├── bond_cleavage_engine.py       #   精确糖苷键切分 + 碳守恒 (15KB)
│   ├── feature_extractor.py          #   Morgan FP + Murcko + 拓扑骨架 (16KB)
│   ├── chemical_classifier.py        #   三层分类引擎 (SMARTS→Glycolipid→Tanimoto) (15KB)
│   ├── molecular_visualizer.py       #   三色高亮 + 分子图渲染 (25KB)
│   ├── stereochemistry_rescue.py     #   ChEMBL/LOTUS/Name rescue (16KB)
│   ├── modification_scanner.py       #   14 SMARTS 修饰扫描 (7KB)
│   ├── secondary_fragment_scanner.py  #   核苷酸/肽键检测 (9KB)
│   ├── taxonomy_lotus_matcher.py      #   LOTUS 物种匹配 (11KB)
│   └── taxonomy_online_resolver.py    #   PubChem/GBIF 在线查询 (20KB)
│
├── scripts/                           # 🔵 执行脚本层
│   ├── run_v13_full_pipeline.py      #   V13 生产管线 (主入口)
│   ├── run_v12_full_pipeline.py      #   V12.3 全量管线 (含 HTML 报告生成)
│   ├── run_glyconp_pipeline.py       #   一键运行入口
│   ├── run_full_refresh.py           #   全量刷新
│   ├── run_benchmark_200_tier_a.py   #   226 条统一基准测试
│   ├── run_audit_fixes.py            #   代码审计修复验证 + CSV 更新
│   ├── rescue_generic_sugars.py      #   NLP 泛指糖救援 Phase 1 (A/C/D)
│   ├── rescue_generic_sugars_phase2.py#  NLP 泛指糖救援 Phase 2 (E/F/G)
│   ├── extract_saponins.py           #   Saponin 子数据库提取 + HTML 报告
│   ├── batch_processor.py            #   批处理框架
│   ├── generate_glyco_landscape.py   #   糖化学图谱可视化 (56KB)
│   ├── generate_saponin_charts.py    #   Saponin 统计图
│   └── nlp_literature_mining.py      #   文献 NLP 挖掘
│
├── tests/                             # 🟡 测试层 (30 files)
│   ├── test_comprehensive_library.py #   参考库完整性测试
│   ├── test_bug_fixes_regression.py  #   回归测试
│   ├── test_precision_cleavage.py    #   精确切分测试
│   └── ...
│
├── data/                              # 数据文件
│   ├── Coconut.csv                   #   COCONUT 全量 (~692MB)
│   ├── benchmark_unified.json        #   226 条统一基准集 (139KB)
│   └── Ginsenoside.csv              #   人参皂苷专题数据
│
├── reports/                           # 输出报告
│   ├── GlycoNP_Deep_Enriched_v13_Final.csv   # V13 终态数据库 (228MB)
│   ├── GlycoNP_Deep_Enriched_v13_Pruned.csv  # V13 修剪版 (163MB)
│   ├── GlycoNP_Saponin_DB_v13.csv           # Saponin 子库 (52MB)
│   └── *.html                                # 可视化报告
│
├── log/                               # 日志
├── docs/                              # 文档/演示
├── CHANGE_LOG.md                      # 变更日志
├── PROJECT_STATUS.md                  # 项目状态报告
└── README.md                          # 本文件
```

---

## 依赖关系规则 (Dependency Rules)

```
✅ scripts/ → lib/     (允许: 脚本调用核心库)
✅ lib/     → lib/     (允许: 库之间互相调用)
⚠️ scripts/ → scripts/ (有限允许: 仅限明确上下游, 如 V13 → V12)
❌ lib/     → scripts/ (禁止: 核心库不得调用脚本层)
```

---

## 依赖 (Dependencies)

| 包 | 版本 | 用途 |
|----|------|------|
| `rdkit` | ≥2023.03 | 分子操作、SMARTS 匹配、指纹、Murcko 骨架 |
| `pandas` | ≥1.5 | 数据表操作 |
| `tqdm` | ≥4.65 | 进度条 |
| `xlsxwriter` | ≥3.0 | Excel 图片嵌入 |
| `numpy` | ≥1.24 | 数值计算 |
| `networkx` | ≥3.0 | 图遍历 (糖链拓扑) |

---

## 快速开始 (Quick Start)

```bash
# 全量运行 V13 管线
python scripts/run_v13_full_pipeline.py

# V12.3 全量管线 (含 HTML 报告, 约 2 小时)
python scripts/run_v12_full_pipeline.py

# 限量测试 (1000 条, 约 20 秒)
python scripts/run_glyconp_pipeline.py --limit 1000

# 226 条统一基准测试
python scripts/run_benchmark_200_tier_a.py

# 代码审计修复验证 + CSV 更新
python scripts/run_audit_fixes.py
```

---

## 基准测试状态 (Benchmark Status)

| 基准集 | 条目数 | 通过率 | 最后验证 |
|--------|--------|--------|----------|
| 统一基准集 (benchmark_unified.json) | 226 | 100% (226/226) | 2026-03-27 |
| Tier A 天然产物 (PubChem) | 50 | 100% (50/50) | 2026-03-25 |
| 1K 分层压力测试 | 900 | 0% 崩溃率, 97.2% 特异性 | 2026-03-26 |

---

## 许可 (License)

本项目仅供学术研究使用。COCONUT 数据遵循其原始许可协议。
