# sun_shallow（任务 3）：Sun 粪便 WGS 浅深度子采样三臂实验

2026-09-12。动机：论文主打 1–2× 测序深度下的 PTR 估计，但真实数据 benchmark 全部
超深度（Sun 每样本 125–135 Gb）。本实验把 3 个 Sun 样本（S01↔SRR13371683、
S06↔SRR13371682、S07↔SRR13371681）子采样到 ~5 Gb 与 ~10 Gb（双端合计，
约为原深度 4% 与 8%），在真实 error + 真实菌群混合下重跑 P2 三臂，
检验 cross-estimator concordance 是否保持、多少物种仍拿到可用估计。

## 1. 子采样（固定种子，确切 read 数已核）

- 工具：`seqtk sample -s42`，_1/_2 两端同一种子同一 N 分别采样（reservoir，
  内存峰值 26–67 GB）；输出 `gzip -1`。
- 目标 N：**5 Gb → 16,666,667 pairs；10 Gb → 33,333,333 pairs**（×300 bp 双端合计）。
- 实测（`subsample_manifest.tsv`，两端 `wc -l` 核实，前 10 万 read 名配对检查
  mismatch=0）：

| srr | depth | n_target | n_1 | n_2 | 原深度占比* |
|---|---|---:|---:|---:|---:|
| SRR13371681 | 5g | 16666667 | 16666667 | 16666667 | 3.69% |
| SRR13371681 | 10g | 33333333 | 33333333 | 33333333 | 7.38% |
| SRR13371682 | 5g | 16666667 | 16666667 | 16666667 | 3.76% |
| SRR13371682 | 10g | 33333333 | 33333333 | 33333333 | 7.52% |
| SRR13371683 | 5g | 16666667 | 16666667 | 16666667 | 3.99% |
| SRR13371683 | 10g | 33333333 | 33333333 | 33333333 | 7.98% |

*原深度按 SRA spots（SRR13371681=451.5M、682=443.5M、683=417.7M）。

- fastq 文件大小（gzip 后）：5g 六端合计 **8.77 GB**（每端 1.44–1.48 GB）；
  10g 六端合计 **17.53 GB**（每端 2.88–2.96 GB）；**总计 ≈26.30 GB**。
  文件保留在 `$B/bench/sun_shallow/fq/`（中间产物不删）。
- 子采样作业（amd，2c/64G）：5g 12–13 min，10g 13–33 min（job 见 §6）。

## 2. 三臂跑法（完全沿用 P2_sigma 口径）

- **arm B sk2bGrow-on-WGS**：`src/target/release/sk2bgrow profile fq1 fq2 -d
  $P2/db_bcgI --mode wms --no-stats --threads 2 --quiet`（BcgI-only 锚库，
  UHGG v2.0 4698 reps，与 Pilea 同一参考集）；stats 用
  `sk2bgrow.cli profile <counts> --db db_bcgI --use-rust-windows --count-model ztp`。
- **arm A Pilea default（C）**：`env-pilea138/bin/pilea profile fq1 fq2 -d
  $P2/pileadb_uhgg -t 16`（自带 gate）。
- **arm A Pilea gates-off（D）**：同上加 `-x 0 -z 0 -c 0`（vendored profile.py
  已打 <2-window 跳过补丁，与既往 D 臂一致）。
- 计时统一 `/usr/bin/time -v`。**注意：既往 P2 全深度跑在 intel 分区，本次按
  现行约定跑在 amd 分区**，跨分区 wall-time 只作量级参考；同深度内部比较同机同线程。
- 评分复用 P2 逻辑（`scripts/ss_50_agreement.py`）：共同分母 = 按 genome merge
  后双侧有限 log2PTR；两档分母——all_reported（与既往一致）与 qc_pass
  （sk2bGrow fraction≥0.5 & dispersion<5）。产物：
  `sun_shallow_agreement.tsv`、`sun_shallow_yield.tsv`。

## 3. yield（拿到估计的物种数）

| sample | depth | sk2bGrow 报告 | sk2bGrow QC-pass | Pilea default | Pilea gates-off |
|---|---|---:|---:|---:|---:|
| S01 | 5g | 329 | 38 | 8 | 680 |
| S01 | 10g | 388 | 63 | 15 | 735 |
| S06 | 5g | 364 | 53 | 14 | 799 |
| S06 | 10g | 465 | 68 | 25 | 832 |
| S07 | 5g | 373 | 68 | 22 | 745 |
| S07 | 10g | 477 | 84 | 32 | 769 |
| S01 | full | 808 | – | 67 | 859 |
| S06 | full | 1019 | – | 80 | 917 |
| S07 | full | 959 | – | 87 | 837 |

- sk2bGrow QC-pass ≈ 全深度的 **8–10%（5g）/ 12–15%（10g）**；Pilea default
  ≈ 全深度的 **12–25%（5g）/ 22–37%（10g）**。真实 1× 深度（~1.3 Gb）会比
  5g 再低一个量级，论文措辞需注意：浅深度下能进共同分母的只有几十个个物种。
- gates-off 产量对深度不敏感（680–832 vs 全深度 837–917）——多出来的"产量"
  依旧是噪声（见 §4）。

## 4. concordance（A vs B，共同分母；x=B log2PTR，y=A log2PTR）

A_default vs B（all_reported）：

| sample | depth | n | r | CCC | bias | LoA |
|---|---|---:|---:|---:|---:|---|
| S01 | 5g | 6 | 0.434 | 0.252 | +0.206 | [−0.19, 0.60] |
| S01 | 10g | 11 | −0.222 | −0.096 | +0.383 | [−0.35, 1.11] |
| S06 | 5g | 14 | 0.144 | 0.105 | +0.153 | [−0.66, 0.96] |
| S06 | 10g | 20 | 0.466 | 0.375 | +0.196 | [−0.38, 0.77] |
| S07 | 5g | 21 | 0.636 | 0.460 | +0.103 | [−0.72, 0.93] |
| S07 | 10g | 28 | 0.564 | 0.429 | +0.134 | [−0.47, 0.74] |
| **ALL** | **5g** | **41** | **0.504** | **0.364** | **+0.136** | **[−0.63, 0.90]** |
| **ALL** | **10g** | **59** | **0.363** | **0.278** | **+0.201** | **[−0.44, 0.84]** |
| S01/S06/S07 | full | 58/68/78 | 0.43/0.34/0.51 | 0.31/0.25/0.37 | +0.25/+0.22/+0.22 | ±~0.4–0.9 |

A_default vs B（qc_pass 分母，ALL）：5g n=38 r=0.420 CCC=0.296；10g n=51
r=0.354 CCC=0.277——与 all_reported 几乎一致（浅的共同分母本来多数就是 QC-pass）。

A_gatesoff vs B：5g ALL n=736 r=−0.047；10g ALL n=834 r=+0.032；qc_pass 分母
r 0.16–0.19——**与全深度结论相同：gates-off 的产量是噪声，Pilea gate 在做真实工作**。

## 5. 结论（回答任务三问）

1. **1–2× 等效深度下 concordance 是否保持？——基本保持。**
   合并 3 样本：5g r=0.50 / CCC=0.36，10g r=0.36 / CCC=0.28，全深度 r=0.34–0.51 /
   CCC=0.25–0.37。系数同量级、bias 略降（+0.14~+0.20 vs 全深度 +0.22~+0.25）、
   LoA 不宽于全深度。即：**深度降到 ~4–8% 后，对仍能同时过两边 gate 的
   几十个个物种，两种独立方法的 log2PTR 一致性没有崩**。
   注意小样本噪声：S01@10g n=11 出现 r=−0.22 的负相关——浅深度单样本共同分母
   只有 6–28 个物种，单样本 r 不可靠，论文只引合并口径。
2. **多少物种仍拿可用估计？** sk2bGrow QC-pass：5g 每样本 38–68、10g 63–84；
   Pilea default：5g 8–22、10g 15–32；两侧都能进的共同分母 5g 6–21、10g 11–28。
   相对全深度（sk2bGrow 报告 808–1019、Pilea 67–87、共同分母 58–78）是
   断崖式下降，且近似随深度线性——符合"PTR 估计需要覆盖度"的基本面，
   也印证 1–2×（~1.3–2.6 Gb）深度下可用物种数会更少。
3. **与全深度对比。** 见 §3/§4 表：concordance 保持、yield 断崖、gates-off
   恒为噪声。对论文的含义：浅深度的瓶颈不是估计一致性，而是**过 gate 的
   物种数**；写 1–2× 场景时应报"每深度可用物种数 + 合并 concordance"两个数。

## 6. 作业与计时记录

- job IDs（`jobs.txt`）：子采样 4025218/4025222/4025226/4025230/4025234/4025238；
  armB 4025219/4025223/4025227/4025231/4025235/4025239；PileaC 4025220/4025224/
  4025228/4025232/4025236/4025240；PileaD 4025221/4025225/4025229/4025233/4025237/
  4025241。24 个作业全部 COMPLETED（2026-09-12 08:27–09:35）。
- armB count（amd，2 线程）：5g 7.7–11.1 min，10g 14.2–19.0 min，RSS ~0.78 GB
  （全深度 intel 6:50）。stats 每样本数分钟内完成。
- Pilea（amd，16 线程）：C 2.2–4.2 min，D 4.0–8.5 min，RSS 0.54–0.74 GB
  （全深度 C 10–77 min、RSS ~1.5–2 GB）。

## 7. 产物

- `fq/`：12 个 subsampled fastq.gz（合计 ≈26.30 GB，保留）。
- `fq/manifest/` + `subsample_manifest.tsv`：种子/确切 read 数/字节数/配对检查。
- `counts_B/<SRR>_<depth>/`：counts.tsv、windows.tsv、stats/output.tsv、time.txt。
- `armA_WGS/{C,D}_<SRR>_<depth>/`：output.tsv、time.txt、pilea.log。
- `sun_shallow_agreement.tsv`、`sun_shallow_yield.tsv`：本文件配套数字。
- 回传本地：scp 至 `data/sun_shallow/`。评分脚本 `scripts/ss_50_agreement.py`，
  重跑：`python scripts/ss_50_agreement.py`（幂等，直接重算两个 tsv）。
