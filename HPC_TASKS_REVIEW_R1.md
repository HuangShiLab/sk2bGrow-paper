# HPC 任务清单 — 模拟评审 R1 回应（P0/P1）

日期：2026-09-11 ｜ 来源： manuscript/reviews/simulated_reviews_2026-09-11.txt
（三位模拟审稿人意见，编号甲=统计、乙=生物学、丙=一致性）

优先级：P0 = 决定论文核心主张，先跑；P1 = 便宜且直接回应最大扣分项。
与既有任务的关系：**σ_eff 已完成；C8 冻结；F7 搁置；screen 口径已定稿（不要重跑 screen）**。
本清单五项全部独立于冻结项。

交付惯例沿用之前：每任务一个 `data/<task_dir>/`（tsv + 数行 REVIEW 结论），
跑完 scp 回论文仓 `data/`，表/图由论文仓本地再生。

---

## 任务 1（P0）F1 对称 harness 重跑 —— 决定 0.5× attribution 是否成立

**评审依据**：甲 M1、丙 M4。同一 arm（anchors+V-fit, k16）同一批 reads，
F1 harness 得 r = 0.568（0.5×，`data/f1_sketch/F1_results.tsv`），
C8 harness（Table 2）得 0.913。两臂 parity flag 不一致
（A 带 `--windows`、E 不带）。0.5× 下"sketch 在存活细胞上反超 panel"
的顺序可能随 harness 翻转——论文最核心的 interaction 主张悬在此上。

**设置**：
- 复用 F1 原 harness 与 reads（Zheng 85-cell 网格，0.5/1/2/5/10×；
  原作业配置见 `data/f1_sketch/REVIEW.md`）。
- 对称化：arm A_k16 与 arm E（sketch, scale 104 与 250）全部统一 flags。
  主跑：两臂都加 `--windows`（对齐当前 harness 中 A 的配置）；
  敏感性：两臂都不带，确认方向不依赖该选择。
- 产物：`data/f1_symm/F1_symm_results.tsv`（列结构同 F1_results.tsv，
  加一列 flag 配置）。

**判据（写进 REVIEW 结论）**：
- A_k16@0.5× 在 F1 harness 的 0.568 是否复现（即差异是 harness 固有还是 flag artifact）；
- 对称条件下 0.5× 的 A-vs-E 排序（存活细胞上谁高）是否翻转；
- ≥1× "两 landmark 源不可区分"的结论是否不动。

**成本**：isolate 级，分钟级/cell，合计 < 2 节点小时。**最优先**——
结果直接决定 §4 的 claim 措辞（interaction vs 融合冗余）。

---

## 任务 2（P0）mm=1 端到端单样本实测 —— 把外推变成测量

**评审依据**：丙 M1。63× lookup 加速与 ~10–25× 全样本收窄
都测于 1M read-pair 子集；审稿人要求至少一个全样本端到端数字。

**设置**：
- 同一 C5 库（24.1M anchor，已在，构建记录见 c5 REVIEW §4d 作业 3995938）。
- count 加 `--max-mismatch 1`（运行时参数，无需重建索引），
  `/usr/bin/time -l` 计时，8 线程，同基线口径。
- 先跑 1 个样本；若 wall 落在外推区间（1.9–3.1 h）再跑满 9 样本数组。
- 基线： `data/c5_review/c5_cost_per_sample.tsv`（16–27 h/样本）。

**产物**：`data/mm1_e2e/mm1_full_sample_cost.tsv`
（样本、wall、峰值 RSS、anchor 保留率、丢 genome 数），
外加 mismatch_hist 对比（确认 5.3% anchor loss 与 no-genome-lost
在全长 reads 上成立）。

**判据**：实测全样本加速比；若显著偏离外推区间，按实测修订
§7c 与 Abstract（当前措辞已标注 extrapolation，改数字即可）。

**成本**：每样本 ~2–3 h wall，9 样本 ~30 节点小时。

---

## 任务 3（P1）Sun 浅深度子采样 —— 在真实数据上检验核心 use case

**评审依据**：乙 M1（生物学意义 2/5 的主因）。论文动机是 1–2× 深度，
但所有真实数据都是超深度（Sun 每样本 125–135 Gb），
动机场景从未在真实 error + 真实菌株混合下检验。

**设置**：
- 源数据：SRR13371682/683/681（2×150，各 125–135 Gb，PRJNA689204，HPC 已有）。
- 子采样到 ~5 Gb 与 ~10 Gb（双端合计，约原深度的 4% 和 8%）：
  `seqtk sample -s <固定种子>`，**记录种子与确切 read 数**。
- 两工具三臂：sk2bGrow（anchors+V-fit）、Pilea defaults、Pilea gates-off；
  参考集与共同分母口径**完全沿用原三臂作业**。
- 每（样本 × 深度）报：共同分母物种数、r/CCC/bias/LoA、
  两工具各自的 yield（拿到估计的物种数与 QC-pass 数）。

**产物**：`data/sun_shallow/sun_shallow_agreement.tsv` + 每深度 yield 表。

**判据**：1–2× 等效深度下 concordance 是否保持、
多少物种仍拿到可用估计——这条直接决定 §7a 能不能从
"concordance validation"升级成"use case demonstrated"。

**成本**：每（样本×深度）约原三臂 run 的 4–8%；合计 < 20 节点小时。

**附项（可选，同一作业顺带）**：Zheng 网格加 1% 替换错误模拟臂
（乙 M1-ii；simulator 已有，只多种一次 reads），闭合
"模拟无测序错误"的 caveat。

---

## 任务 4（P1）dedup 机制回归 —— 把"不要 dedup"从经验变成机制

**评审依据**：乙 M3。dedup 推荐是湿实验用户最可能引用的一条，
但目前没有机制解释，也没有超出"一个库一个酶一个中心"的证据。

**设置**：
- 输入：Sun 真实 2bRAD 库（BcgI）per-anchor 原始计数
  （HPC 上 sun_three_arm 管线的 count 输出，未 dedup 版）。
- 计算 per-anchor duplication fraction（raw count / exact-dedup count，
  或以 UMI 等价标记），对每个物种回归于 circular distance from ori。
- 机制假设：PCR 偏好高拷贝模板 → ori 近端（复制中多拷贝）duplicate
  fraction 更高；exact dedup 恰好除掉的正是梯度本身。
- 样本 n = 3，别过度检验：报每物种斜率、符号一致性、 anchors n。

**产物**：`data/dedup_mech/dedup_mechanism.tsv`（物种 × 斜率 × n × R²）
+ 三行机制结论 + 边界条件声明（cycle 数/input mass 未知；对其他 2bRAD
protocol 的可迁移性是 future work）。

**判据**：3/3 样本斜率为正即支持机制。

**成本**：纯分析，< 1 节点小时。

---

## 任务 5（P1）bootstrap CI + 2×2 交互形式检验 —— 本地，不占机时

**评审依据**：甲 M3/M9。全文无一处不确定度量化；所有"within noise"
"indistinguishable"都是目测；2×2 交互从未作为交互被检验。

**设置（在论文仓本地跑，不需要 HPC）**：
- 输入已在仓内：`data/results_raw.tsv`（Table 2）、
  `data/f1_sketch/F1_results.tsv`（F1）、`data/panel_sweep.tsv`（Table 7）。
- media bootstrap（16 media 重抽 10k 次）给 r/slope/RMSE 的 95% CI；
- Fisher z 配对检验（同一 harness 内 A vs E，逐深度）；
- 产出每张表加 CI 列所需的中间文件 `data/ci_bootstrap/`。

**判据**："≥1× 两 landmark 源不可区分"量化为 CI 重叠；
"r 在 k=8 达峰"给出 CI 后是否仍成立；为 Methods 的
pre-specified vs post-hoc 声明提供素材。

**成本**：纯重分析，笔记本分钟级。

---

## 顺序与依赖

1. **任务 1 最优先**（可能改变 §4 claim 与多个数字，先拿结果再定稿）。
2. 任务 2/3/4 互相独立，可并行提交。
3. 任务 5 本地随时可做（我可以直接在论文仓跑）。
4. 全部完成后：Table 2/7 加 CI 列、§4/§7c 数字按新结果修订、
   图 6/10/12 必要时重绘。

## 口径红线（本清单各任务通用）

- 禁跨 harness 比较数字（任务 1 的存在就是为了消灭这种比较）。
- 子采样必须固定种子并记录确切 read 数。
- 计时比较必须同机、同线程数、`/usr/bin/time -l` 同口径。
- Sun/C5 报数一律共同分母 + QC-pass；禁报裸 recall。
- 不重跑 screen、C8、σ_eff——这些已定稿。

## P2 暂缓（记录备查，本次不跑）

- F2 加低 GC、多 Mb 基因组解基因组大小混杂（甲 M6/乙 M5）
- Pilea full-depth 直跑对齐 0.976（乙 M6；~55 GB）
- 45k *Escherichia* assembly 的 ANI×N50 子集 sweep（丙 top-3）
- C5 用现行 stats 层重融合（QC 列审计，丙 M7）
- 粪便物种 log2PTR 排序 vs 生理学先验的对比分析（乙 M2-i，纯重分析）
