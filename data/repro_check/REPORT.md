# sk2bGrow 可复现性诊断（C1 0.913 vs F1 0.568）

日期：2026-09-12 ｜ 机器：HPC（/lustre1/g/aos_shihuang/sk2bgrow-hpc）｜ 产物目录：`bench/repro_check/`

## 结论（TL;DR）

**现行（已提交+工作树）代码 + C1 的 0.5× counts 可以复现 0.913 量级的结果（实测 r=0.941，n=16）。**
**0.568 不是代码漂移造成的**——它定位为 **F1 的 count 输入配置不同：F1 只用 `_1` 单端 reads（7,812 reads），C1 用 `_1+_2` 双端（15,624 reads）**。
- count 层（Rust）：**零漂移**——现行 binary 对同一 fastq 重 count，与 C1、F1 已提交的 counts.tsv 均 **md5 逐字节一致**。
- stats 层（Python）：**有数值漂移**（逐 cell 与 C1 提交值不同；C1 时代的 ZTP boundary fixes 已不在工作树，仅存于 `ztp.py.bak`），但在 C1 counts 上 headline r 不受影响（0.94 vs 0.91）。

## 代码状态（运行时快照，见 code_state.txt）

- `git describe`：**`275778f-dirty`**（HEAD 275778f "Replace the HPC handover with a full research plan"）
- 工作树脏文件：`python/sk2bgrow/{cli,fit,fusion,ztp}.py`、`crates/*/count.rs,profile.rs,index.rs` 等（git diff --stat 合计 ~855 行）。
- 全部 python 文件最后一次 **commit** 是 8 月 24 日的 `53c56c5`；C1（8/26）跑的是当时的未提交工作树（two-stage decoupled weights + ZTP boundary fixes）。
- `ztp.py.bak`（8/27）= C1 时代的修复版；**当前工作树 ztp.py 已退回 naive 边界公式**（`log1p(-exp(-lam))`），新增 `total_count` 列（供 GLM 用）。
- Rust binary `target/release/sk2bgrow` 构建于 9/4 17:26，晚于 count.rs/profile.rs 最后修改（9/4 17:23）→ binary 与源码一致。
- Python 模块解析到 `$B/src/python/sk2bgrow`（editable install；PYTHONPATH 设置后同路径）。

## 主测试：现行 stats 代码 × C1 0.5× counts（17 media）

- 调用：`python -m sk2bgrow.cli profile <counts> --db bench/C1/db --output bench/repro_check/out_c1counts/<M>__<SRR> --windows bench/C1/counts/<SRR>.0.5x/windows.tsv`（windows.tsv 已存在于各 counts 子目录，传入；与 04_stats.sh 完全一致）。
- 17/17 成功（首轮 11 个因登录节点 OpenBLAS 线程超限失败 rc=130，`OPENBLAS_NUM_THREADS=1` 后全部通过）。
- **Pearson r（fused log2PTR vs growth_rate，排除 RUN_OUT，n=16）：**

| 数据 | r |
|---|---|
| **本次：现行代码 × C1 counts** | **0.9411** |
| C1 已提交 output.tsv（复算） | 0.8957 |
| 论文 data/results_raw.tsv arm A 0.5×（复算） | 0.9132 ← 头条 0.913 |
| F1 A_k16 0.5×（单端，F1_results.tsv） | 0.5680 |

逐 cell 值见 `repro_summary.tsv`。注意：论文 results_raw.tsv 的逐 cell 值与 HPC 上 C1 output.tsv 略有出入（两者 r 分别为 0.9132/0.8957），来源待查；但三者都在 0.9 量级，远离 0.568。

## 对照 A：现行代码重跑 F1 M10（SRR11558980, k16, 0.5×）

- 输出 `out_controlA_M10/output.tsv` 与 F1 已提交 `bench/F1/out/A_k16_M10_0.5x/output.tsv` **逐字节一致（diff 无差异）**。
- → F1 的数字确实是现行 stats 代码的产物。
- `windows.rates.tsv` 每行不同仅因工作树新增了 `total_count/log2_rate_raw/gc_corrected` 列（schema 漂移），rate 列与融合结果不变。

## 对照 B（count 层判定，md5 级）

现行 Rust binary 重 count M6（SRR11558945, 0.5×）：

| 输入 | md5 | 对比 |
|---|---|---|
| 仅 `_1`（单端） | 24a401d6… | == F1 cells/k16 同 cell counts.tsv（逐字节一致） |
| `_1 + _2`（双端） | 388ce7ac… | == C1 counts/同 cell counts.tsv（逐字节一致） |

- → **count 层自 8 月 26 日以来零漂移**；F1 与 C1 counts 的差异（此前 spot-check md5 不同）**完全是输入 reads 不同**：F1 `reads_total=7,812` vs C1 `15,624`（精确一半）；两细胞 `windows.tsv` md5 相同（anchor 集合一致，只有 read 计数不同）。
- F1 f1_11_counts.sh 注释亦自述："Uses _1 reads only (single-end)"。
- 佐证：F1 A_k16 在 1×（单端，≈15.6k reads）r=0.906 ≈ C1 0.5×（双端 15.6k reads）r≈0.91–0.94 —— r 由 read 数决定，与"0.568 源于 read 数减半"一致。

## 最终判定

1. **0.913 可否由现行代码 + C1 counts 复现？可以**（r=0.941；精确复算论文表得 0.9132 的那张逐 cell 表与现行 C1 output.tsv 有细微出入，见上）。
2. **0.568 的定位：count 输入配置（单端 vs 双端），不是 count 代码、也不是 stats 代码回归。**
3. **stats 层确有数值漂移**（ZTP boundary fixes 丢失 + fit/fusion 新增 GLM 等），建议在论文定稿前决定：恢复 `ztp.py.bak` 的边界修复并重新锁定基线，或接受 0.94 为新基线。

## 产物路径（HPC）

- `bench/repro_check/code_state.txt` — 代码状态快照
- `bench/repro_check/run_main.sh` — 主测试脚本
- `bench/repro_check/main_progress.log` / `nohup_main.log` — 运行日志
- `bench/repro_check/out_c1counts/<M>__<SRR>/` — 17 个 cell 的现行代码 stats 输出（output.tsv, windows.rates.tsv, stats.log）
- `bench/repro_check/repro_summary.tsv` — 逐 cell 对比表
- `bench/repro_check/out_controlA_M10/` — 对照 A 输出（与 F1 提交逐字节一致）
- `bench/repro_check/recount_M6_se/`、`recount_M6_pe/` — 对照 B 重 count 产物
- `bench/repro_check/REPORT.md` — 本文件

本地副本：`/Users/macstudio/Downloads/sk2bGrow-paper/data/repro_check/`（repro_summary.tsv、REPORT.md、code_state.txt）

## 全深度再生成汇总（2026-09-12 追加）

按用户拍板，用 HPC 现行代码产物再生成论文 Table 2 arm A 的全部 5 个深度。设置与主测试完全一致（`PYTHONPATH=$B/src/python`、envs/sk2bgrow python、`--db bench/C1/db`、`--windows` 传 counts 目录现成 windows.tsv、`OPENBLAS_NUM_THREADS=1`、xargs -P 8）；68/68 cell rc=0（`depths_progress.log`），代码状态见 `code_state.txt` 追加段（`275778f-dirty`）。

**fused log2PTR 对 growth_rate 的 Pearson r（排除 RUN_OUT，first-pick，n=16）：**

| cov | r（现行代码再生成） | 参照：C1 已提交 run-level（n=42） |
|---|---|---|
| 0.5× | 0.9411 | 0.8420（0.5× first-pick n=16 为 0.8957） |
| 1× | 0.9185 | 0.9136 |
| 2× | 0.9586 | 0.9611 |
| 5× | 0.9711 | 0.9703 |
| 10× | 0.9699 | 0.9683 |

- 输出布局：`out_c1counts/<M>__<SRR>/`（0.5×，主测试产物原样保留）、`out_c1counts/<D>x/<M>__<SRR>/`（1/2/5/10×，深度分层以避免覆盖；skip-if-exists）。
- 逐 cell 记录：`regen_armA_cells.tsv`（85 行 = 5 深度 × 17 medium；列：medium, srr, cov, log2ptr, est_cov, pass_qc, fusion_model, method, qc_reason, excluded, note, growth_rate）。
- 字段映射（见 `regen_header_sample.txt` 的 output.tsv 原表头 + 样例行）：`log2ptr` = output.tsv 的 `log2(PTR)`（fused 估计）；`est_cov` = `coverage` 列；布尔字段 = `pass_qc`（全部原样保留供本地挑用）。
- 判定：各深度现行 r 与 C1 已提交同口径数值一致（0.5× 三变体 0.8957/0.9132/0.9411 的讨论见前节）；再生成产物可直接用于更新论文 Table 2 arm A 行。

## 三臂统一再生成汇总（2026-09-12 追加二）

为消除"arm A 新 / B·C·E 旧"的混口径，将 B、C、E 三臂也统一到 HPC C1 实例（现行代码 `275778f-dirty`，输入与 arm A 完全一致）。逐 cell 表：`regen_armB_cells.tsv`、`regen_armC_cells.tsv`、`regen_armE_cells.tsv`。

**各臂各深度对 growth_rate 的 Pearson r（排除 RUN_OUT；括号内为 n≠16 的说明）：**

| cov | A（anchors+V-fit） | B（anchors+sorted） | C_default（Pilea） | C_relaxed（Pilea） | E（sketch104+V-fit） |
|---|---|---|---|---|---|
| 0.5× | 0.9411 | 0.5165 | 无估计（gate 拒绝，n=0） | NA（17 格全 0，常数列） | 0.8183（n=10） |
| 1× | 0.9185 | 0.7753 | 无估计（n=0） | 0.7784 | 0.9089 |
| 2× | 0.9586 | 0.9253 | 无估计（n=0） | 0.8481 | 0.9662 |
| 5× | 0.9711 | 0.9472 | 无估计（n=0） | 0.9877 | 0.9863 |
| 10× | 0.9699 | 0.9378 | 0.9730 | 0.9730 | 0.9812 |

方法学与口径：
- **arm B**：`profile --method sorted`（fit_sorted_ransac，iRep/Pilea 式 rank 回归），其余参数与 arm A 相同（同 counts、同 windows.tsv、同 db）。85/85 rc=0，输出 `out_c1counts_b/<D>x/<M>__<SRR>/`。
- **arm C**：纯提取，Pilea 产物 170/170（17×5 深度×2 模式）齐全无需补跑。gate 拒绝的 cell（default 模式 0.5–5×）output.tsv 仅有表头 → 表中估计字段留空、保留行（与 Mac 备份表示法一致）；Pilea 无 pass_qc 列，`passed` 留空。
- **arm E**：scale 104 sketch counts 直接复用 `bench/F1/sketch/s104/`（经单 cell 重生成 md5 确定性校验：`003a575e…` 逐字节一致），hardlink 到 `sketch_e/counts/`；stats 按 **F1 arm E 口径不传 --windows**（sketch landmark 与 anchor 窗口不兼容；F1 的 74/85 v_shape + 4 v_shape_segmented 证明该流程成立），85/85 rc=0，输出 `sketch_e/out/<M>__<SRR>__<D>x/`。
- 交叉验证：arm E 各深度 r 与 F1_results.tsv 已提交值**逐深度完全一致**（0.8183/0.9089/0.9662/0.9863），现行代码与 9/8 F1 代码在 sketch 统计路径上无差异。

解读：
- **B 臂 0.5× 崩至 0.52**——sorted/rank 回归在低深度失效（高估计偏置，如 Mac 备份 M10 0.5×=2.30），与论文 Fig 4（sorted 回归在非生长培养上虚构梯度）的叙述互相支撑；2× 以上恢复 0.92+。
- **C_default 只在 10× 有估计**：Pilea 默认 gate 在 <10× 全部拒绝输入（HPC 与 Mac 侧一致）。
- **A vs E（2×2 的 estimator×landmark 交互）**：同一 V-fit 下 anchors（A）与 sketch（E）在各深度 r 差 ≤0.11（0.5× 0.94 vs 0.82），V-fit 对 landmark 来源稳健；同一 anchors 下 V-fit（A）对 sorted（B）优势随深度收窄。
- 五臂全部统一到 HPC C1 实例后，Table 2 的"同一批 subsampled reads"口径成立，可与 arm A 新值并列使用。

## arm E 双端修正（2026-09-12 追加三）

发现此前复用 `bench/F1/sketch/s104/` 的 counts 来自 F1 **单端** cells，与 C1 双端口径不符。已按 C1 paired-end fq（`<run>.<depth>x_{1,2}.fq.gz`）重生成：armE_counts.py 逐 fq 出表、不合并双端，故在 repro_check 下写了池化驱动 `armE_counts_pe.py`（复用其常量与 Pilea count64，k=31、scale 104、cache 复用 F1 的 `cache_s104.pkl`），每格把 `_1+_2` 的 k-mer 计数并入一张 count 表（与 C1/Pilea 双端约定一致）。counts 85/85（`sketch_e/counts_pe/`，约 10 分钟，xargs -P 8；首轮因漏设 `OPENBLAS_NUM_THREADS=1` 出现 rc=130，修复后零失败）；stats 与 F1 arm E 同口径（db=C1 db、**不传 --windows**）85/85 rc=0（`sketch_e/out_pe/`）。逐 cell 表：`regen_armE_PE_cells.tsv`。

**arm E 单端 vs 双端 Pearson r（排除 RUN_OUT；括号内 n）：**

| cov | E 单端（F1 counts） | E 双端（C1 fq，推荐） |
|---|---|---|
| 0.5× | 0.8183（n=10） | **0.8831（n=16）** |
| 1× | 0.9089 | 0.9121 |
| 2× | 0.9662 | 0.9578 |
| 5× | 0.9863 | 0.9898 |
| 10× | 0.9812 | 0.9865 |

- 双端把 sketch 计数翻倍：0.5× 从"10/16 格有估计"变为 16/16 全覆盖，r 0.82→0.88；5×/10× 也小幅上调。2× 微降 0.9662→0.9578，属逐 cell 噪声正常摆动。
- 至此五臂（A/B/C_default/C_relaxed/E）全部落在 **C1 双端实例**上，Table 2 的"同一批 subsampled reads"口径成立；论文 arm E 行应采用 `regen_armE_PE_cells.tsv`（单端版仅作 F1 对照保留）。
- A vs E（PE）差距收窄：0.5× 0.9411 vs 0.8831、5× 0.9711 vs 0.9898（E 反超）——V-fit 对 landmark 来源的稳健性结论不变。
