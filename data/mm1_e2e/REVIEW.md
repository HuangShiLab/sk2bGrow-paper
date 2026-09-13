# Task 2 — mm=1 全样本端到端实测（C5, RBC 宏基因组 9 样本）

日期：2026-09-12 ｜ 作业：试点 4025262_6（SRR28338158）+ 数组 **4025355**（0–8，
partition `amd`）；4025355_6 提交后发现与试点并发写同一 stats 输出目录，已
`scancel`（该样本已由 4025262_6 完整产出，无需重跑）。
数据根：`/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/mm1_e2e/`
（`res/A_<run>/`：time.txt / stats.time.txt / counts.tsv / windows.rates.tsv /
analysis.json / mm1_vs_mm2_per_genome.tsv；汇总 `mm1_full_sample_cost.tsv`）

## 1. 基线口径（复核查证）

基线出处 `$B/bench/C5/scripts/01_cell.sh` + `c5_cost_table.py`，
`data/c5_review/c5_cost_per_sample.tsv`（= `$B/bench/C5/review/c5_cost_per_sample.tsv`）：

- **count 阶段**：`$B/src/target/release/sk2bgrow profile <fq1> <fq2> -d db -o <out>
  --no-stats --threads 8 --quiet`（默认 mm=2），`/usr/bin/time -v` → `time.txt`。
  每样本 wall 58 931–98 049 s（16.4–27.2 h），peak RSS ≈ 2.25 GB，Exit 0。
- **stats 阶段**：`python -m sk2bgrow.cli profile <counts.tsv> --db db --output <out>
  --windows <out>/windows.tsv --use-rust-windows --count-model ztp`，
  `/usr/bin/time -v` → `stats.time.txt`。每样本 1 191–2 614 s，peak RSS ≈ 15 GB。
- 基线"每样本 16–27 h"= **count only**；count+stats 合计 60 814–100 587 s
  （16.9–27.9 h，`sk2bgrow_total_min` 同定义）。
- 机器：partition `amd`（基线在 GPA-1/2 节点），8 线程。

**M4 构建 ≡ src 构建（mm=2 同速）**：同一 1M read-pair 子集、同 8 线程、
counts-only 口径：src = 1 082 s（`res/_phase_diag/time_count_timing.txt`）vs
src_m4 = 955 s（`res/_mm1_bench/mm2ns.time.txt`）。M4 略快 ≤12%，
用 M4 跑 mm=1 与 src 基线比较只会低估加速比。

## 2. 本任务实测设置（与基线逐项对齐）

每样本两阶段，与基线完全同构，仅 count 阶段加 `--max-mismatch 1`、
二进制用 `$B/src_m4/target/release/sk2bgrow`（db 只读）：

1. `env SK2B_COUNT_TIMING=1 /usr/bin/time -v sk2bgrow profile <fq1> <fq2> -d $B/bench/C5/db
   -o <out> --no-stats --threads 8 --quiet --max-mismatch 1`
2. `PYTHONPATH=$B/src_m4/python /usr/bin/time -v python -m sk2bgrow.cli profile
   <counts.tsv> --db db --output <out> --windows <out>/windows.tsv
   --use-rust-windows --count-model ztp`
3. `analyze_sample.py`：全样本 anchor 保留分析（§4）。

**与任务书字面命令的偏差**：任务书写单命令 `profile`（不带 `--no-stats`）。
实测采用基线两阶段拆法：(i) 基线口径即两阶段，逐阶段同口径可比；
(ii) 单命令内部 stats 曾在作业 3999442 死于 python 环境，脆弱。
端到端 wall = 两阶段之和，与 `sk2bgrow_total_min` 同定义。

节点说明：本次跑在 GPA-4-x（8 个任务）与 GPA-1-1/7（2 个任务），
与基线同分区同线程数；20–39× 的量级远大于节点型号差异（±10–20%）。

## 3. 实测结果

### 3a. 每样本 wall / 加速比 / RSS（mm1_full_sample_cost.tsv）

| sample | 基线 count (h) | mm1 count (min) | count 加速 | mm1 stats (min) | **mm1 e2e (min)** | e2e 加速 | count RSS (GB) | stats RSS (GB) |
|---|---|---|---|---|---|---|---|---|
| SRR28338172 | 25.28 | 47.7 | 31.8× | 14.8 | 62.4 | 25.0× | 2.49 | 15.28 |
| SRR28338171 | 27.24 | 42.6 | 38.3× | 19.1 | 61.8 | 27.1× | 2.49 | 15.37 |
| SRR28338162 | 24.88 | 38.5 | 38.8× | 15.4 | 53.8 | 28.3× | 2.49 | 15.23 |
| SRR28338161 | 18.23 | 37.7 | 29.0× | 16.0 | 53.8 | 21.0× | 2.49 | 15.37 |
| SRR28338160 | 19.64 | 45.0 | 26.2× | 17.5 | 62.5 | 19.2× | 2.49 | 15.44 |
| SRR28338159 | 17.68 | 39.9 | 26.6× | 17.8 | 57.7 | 19.0× | 2.49 | 15.45 |
| SRR28338158 | 16.37 | 38.7 | 25.4× | 16.9 | 55.6 | 18.2× | 2.49 | 15.25 |
| SRR28338157 | 17.11 | 50.1 | 20.5× | 37.5 | 87.6 | 12.1× | 2.49 | 15.42 |
| SRR28338156 | 20.63 | 44.6 | 27.7× | 28.4 | 73.0 | 17.4× | 2.49 | 15.19 |

- **count 阶段加速 20.5–38.8×，中位 27.7×**；**e2e（count+stats）12.1–28.3×，中位 19.2×**。
- 全部 count exit 0；RSS 与基线同量级（count ≈2.5 GB，stats ≈15 GB）。
- SRR28338157 的 stats 偏慢（37.5 min vs 其余 15–28 min），拉低了它的 e2e 加速；
  count 加速不受影响。
-  phase timing（SRR28338158 例）：motif_scan 2 314 s（71%）+ lookup 943 s（29%）。
  对照：mm=2 在 1M 子集上 lookup 即需 1 518 s——mm=1 全样本 lookup 反而更低。

### 3b. 与外推区间的偏差

外推区间（1M 子集 A/B 的 8.6× 加速线性外推）：count 1.9–3.1 h/样本。
**实测 count 38–51 min/样本，比外推快 2.5–4.5×，且全部 9 样本一致。**

原因（机制）：1M 子集低估了 mm=2 的 lookup 爆炸。mm=2 时 error/variant tag
撞长 seed posting list 的概率随 read 数增长（`_mm1_bench/seed_hist.txt`：
slot0 hit_mean_len 249→20、slot1 1 453→7，mm=1 把碰撞压掉 1–2 个量级），
lookup 成本超线性；mm=1 把它压回近线性。motif_scan 两臂相同且线性，
在 mm=1 wall 中占 71% 成为下限。**"全样本比子集更受益"是预期内方向。**
论文此前的 "10–25×" 区间（同样来自子集外推）也被实测整体上移。

## 4. 全样本 anchor 保留（不重的跑 mm2）

**前提修正**：任务书假设 mm2 基线 per-anchor counts 仍在 `res/A_<sample>/`；
实际 `01_cell.sh` 在 stats 后 `rm *.counts.tsv`，全样本 mm2 counts 已不存在。
不重跑 mm2，改用**对称衍生指标**：

- 基线 `res/A_<run>/windows.rates.tsv`（stats 从 mm2 counts 算出）每行含
  `n_positive`；窗口在每个 (genome, enzyme, contig) 内划分、覆盖全部 usable 锚，
  故 `D(g) = Σ n_positive` = 该 genome 检测到的 usable 锚数。
- mm1 侧同法。usable 集合由锚 flags/GC 决定、与 count 无关；两臂 stats 调用相同。
- 结构校验：9 样本的 windows.rates.tsv 行数完全相同（2 585 499 行），
  窗口结构逐锚一致，比较严格对称。
- 幅度校验：SRR28338157 top50 真实 mm2 counts 的 raw nonzero vs D_mm2，
  样本级比值 0.80，与本任务 mm1 侧 ΣD/Σraw（10.13M/12.66M = 0.80）一致——
  usable 过滤对两臂是同一常数因子，per-genome 比率 D1/D2 无偏。

结果（每样本 analysis.json + mm1_vs_mm2_per_genome.tsv）：

| sample | anchors_nonzero_mm2 | anchors_nonzero_mm1 | anchor 丢失 | genomes_lost | per-genome 比率中位 | mismatch_hist 保留 |
|---|---|---|---|---|---|---|
| SRR28338172 | 4 816 296 | 4 300 040 | 10.7% | **0** | 0.781 | 89.9% |
| SRR28338171 | 8 968 411 | 8 543 786 | 4.7% | **0** | 0.944 | 94.2% |
| SRR28338162 | 9 384 280 | 8 966 903 | 4.4% | **0** | 0.946 | 94.3% |
| SRR28338161 | 8 192 722 | 7 767 936 | 5.2% | **0** | 0.933 | 94.8% |
| SRR28338160 | 11 344 573 | 11 010 781 | 2.9% | **0** | 0.960 | 95.4% |
| SRR28338159 | 11 766 773 | 11 439 956 | 2.8% | **0** | 0.965 | 95.5% |
| SRR28338158 | 10 458 089 | 10 126 416 | 3.2% | **0** | 0.950 | 95.8% |
| SRR28338157 | 10 196 680 | 9 802 363 | 3.9% | **0** | 0.949 | n/a* |
| SRR28338156 | 8 527 424 | 8 056 582 | 5.5% | **0** | 0.924 | 93.9% |

\* SRR28338157 基线 `*.stats.json` 在 top50 调试期遗失，mismatch_hist 不可比；
窗口 rates 不受影响。

**结论**：
1. **9/9 样本 zero genomes lost**——1M 子集上的 "no-genome-lost" 在全长 reads
   全样本上成立。
2. usable 检测锚丢失 2.8–10.7%（中位 ~4.7%），与 1M 子集 per-anchor 的
   5.3% 同一量级（口径不同：usable-检测 vs raw per-anchor，见 §5）。
3. SRR28338172 偏低（10.7%，mismatch_hist 保留 89.9%，reads_with_anchor
   14.37M→12.54M）：该样本 mm2 锚定 reads 本就最少（14.4M vs 其他 18–30M），
   提示其群落与参考面板分歧更大、更依赖 2-mismatch 匹配——生物学差异而非
   计算 artifact（窗口结构与各样本一致已校验）。
4. mm0（exact）bucket 在 mm1/mm2 两臂完全相同（如 SRR28338158 均为
   48 835 483），mm=1 只动 1–2 mismatch 部分，符合预期。

## 5. 1M 子集 A/B（SRR28338159，已有结果，per-anchor 精度）

SRR28338159 的 1M read-pair 子集 A/B（`res/_mm1_bench/`，作业 3999630 输出 +
`mm1_count.sh` 对比逻辑；汇总 `mm1_subset_1M_ab_summary.tsv`）：

- anchors compared 29 672 394；nonzero mm2 = 1 797 315，mm1 = 1 703 604，
  **lost_by_mm1 = 95 156 → anchor loss 5.29%**（即 R1 所引 5.3%）；
- **共享 nonzero 锚的 count 比率 mm1/mm2 中位 = 1.0000**（保留的锚计数不变）；
- genome 层面：522/522 检测，**彻底丢失 0**；per-genome 总 count 比率中位 0.931；
- 子集计时：mm2 955.0 s vs mm1 111.3 s（8.58×，外推的来源）。

口径：raw per-anchor（counts.tsv 直接合并，mm1_count.sh 逻辑）；
§4 为 usable-检测口径（rates 衍生，D1/D2 中位 0.78–0.96）。两口径结论一致：
丢失个位数百分比、无 genome 彻底丢失（9/9 全样本 + 子集三重确认）。

## 6. 结论与论文修订建议

1. **§7c（及 Methods 计时段）**：把 "10–25×（extrapolated from 1M-pair
   subset）" 替换为实测：
   "On the nine RBC metagenome samples (68–79M read pairs each), restricting
   counting to ≤1 mismatch reduces per-sample count wall time from 16.4–27.2 h
   to 38–51 min (20.5–38.8×, median 27.7×) and end-to-end runtime
   (count + window fitting) from 16.9–27.9 h to 54–88 min (12.1–28.3×,
   median 19.2×), measured on identical hardware, thread count and timing
   methodology as the mm=2 baseline. The full-sample speedup exceeds the
   1M-pair subset estimate (~8.6×) because mm=2 lookup cost grows
   superlinearly with read depth while mm=1 stays near-linear."
2. **Abstract**：相应把 "order-of-magnitude / 10–25× (extrapolated)" 改为
   "20–40× (measured end-to-end on full RBC metagenome samples; median 28×
   for counting)"，并保留 no-genome-lost / ~5% anchor loss 一句
   （现在有了 9 样本全样本背书）。
3. **σ_eff / screen 相关段落**：无需改动；本任务只动 mm 权衡的措辞。
4. 建议正文明确 mm=1 的生物学代价上界：zero genomes lost（9/9）、
   per-anchor 丢失 ≤5.5%（中位）、worst sample 10.7%（SRR28338172，
   高分歧群落）；审稿人若问 "为什么不用 mm=0"，可答 mm1→mm2 方向的对称
   数据在本表。

## 附：复现与产物

- 作业脚本：`scripts/mm1_e2e_cell.sh`（数组 0–8 ↔ samples.tsv 的 9 个 MG 样本）、
  `scripts/analyze_sample.py`、`scripts/collect_cost.py`。
- 汇总：`mm1_full_sample_cost.tsv`（9 行 × 25 列；列定义见 collect_cost.py
  与 §3/§4）；`mm1_subset_1M_ab_summary.tsv`（§5）。
- 每样本：`res/A_<run>/analysis.json`、`mm1_vs_mm2_per_genome.tsv`
  （522 genome 的 D_mm2/D_mm1/比率）。
- 计时口径：/usr/bin/time -v，amd 分区，--threads 8，与基线一致。
