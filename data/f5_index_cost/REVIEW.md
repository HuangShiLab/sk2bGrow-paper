# F5 REVIEW — 匹配密度下两种模式的索引大小与成本（喂 C6/GTDB 预算）

日期：2026-09-09。任务：[`HPC_TASKS_SKETCH_MODE.md`](../../src/benches/HPC_TASKS_SKETCH_MODE.md) F5。
数据：`F5_index_cost.tsv`（1514 行逐条测量）、`F5_half_density.tsv`（60 cell 聚合）、
`f5_scales.tsv`（300 基因组匹配 scale 表）、`genomes.tsv`。产物：`genomes/`（300）、
`db/`（600 单基因组索引）、`db300/`（两个 300 基因组合并索引）、`sketch/`（900 个
FMH 参考 sketch）、`fq/`、`cells*/`、`out/`、`fq_real/`、`cells_real/`、`scripts/`。

作业链（intel partition，全部 ExitCode 0）：select 4010656（2m03s）→ index 4010664
（6m29s）→ sketch 4010703（2m43s）→ 并行 profile_cost 4010764（20m19s）与
half_density 4010732+4010778（22m03s + 11m19s 补 sketch 臂 stats）。中间失败
4010650/4010681/4010694/4010731/4010745 均为脚本 bug（明文 FASTA、导出变量 typo、
`set -o pipefail` 与 zcat|head 冲突、sk2bgrow profile 默认拉起统计层），已修复，
产物幂等无污染。

## 0. 口径（必须随结果一起报）

1. **基因组**：本地 GTDB R226 池（`/lustre1/g/aos_shihuang/databases/GTDB/GTDBr226`，
   732,309 组装）随机抽样 2200、按 GC 均匀取 282，加 F2 的 18，共 **300 个**
   （1–20 Mb，GC 24.3–73.2%，实测自 FASTA），stage 到 `bench/F5/genomes/`。
   合并索引 `db300_*` 一次喂入 300 个 FASTA。
2. **酶臂**：`sk2bgrow index`（release binary），k16 = 全酶板、k8 = F1 top-8。
   **FMH 臂**：binary 尚无 `--mode fracminhash`，sketch 臂的"索引"按其现行实现
   计 = Pilea 式参考 sketch（`sketch_loci`，armE 同款哈希语义，k=31），
   scale 匹配沿用 F2 口径 `scale = round(L / n_anchors)`；s200 为半密度臂。
3. **测量**：构建均包 `/usr/bin/time -v`（wall + 峰值 RSS）；磁盘 = `du -sb` /
   pkl 字节；profile 用 PRJNA1280254 两个真实 WGS 样本各取前 2M 条单端 reads
   （8 线程，`--no-stats`，与 F2 profile 成本口径一致）。
4. **300 基因组合计**：1.048 Gbp，k16 锚点 9,979,215（**9,519/Mb**），
   k8 7,546/Mb；匹配 scale 中位 107（65–262）——**GTDB 平均口径下
   "scale 200 ≈ k16 密度的一半"成立**（实测 s200/k16-匹配 landmark 比
   中位 0.534，范围 0.33–1.33；>1 者均为低 GC 基因组，其匹配 scale 本就 >200）。

## 1. 索引成本对照（F5_index_cost.tsv）

**单基因组构建（300 个的中位）**：

| 模式 | build wall | 峰值 RSS | 磁盘 B/landmark | landmarks 中位 |
|---|---:|---:|---:|---:|
| 酶板 k16 | 3.0 s | 11.2 MB | 38.7 | 27,692 |
| 酶板 k8 | 1.0 s | 10.1 MB | 38.9 | 22,150 |
| FMH 匹配 k16（pkl） | ~1–2 s | ~140 MB | 23.4 | ~30,000 |
| FMH s200 半密度（pkl） | 1.0 s | 137 MB | 23.3 | 14,586 |

**合并索引 db300（GTDB 式一次建库）**：

| 库 | anchors | build wall | 峰值 RSS | 磁盘 |
|---|---:|---:|---:|---:|
| db300_k16 | 9.98 M | **86 s** | **3.93 GB** | **385.9 MB**（38.7 B/anchor） |
| db300_k8 | 7.91 M | 94 s | 3.50 GB | 307.3 MB（38.9 B/anchor） |

**单样本 profile（2M 真实 reads，8 线程）**：

| 对象 | wall | 峰值 RSS |
|---|---:|---:|
| db300_k16 | 530 s / 231 s（两个样本） | ~920 MB |
| db300_k8 | 155 s / 140 s | ~755 MB |
| 单基因组 k16（g08） | 74 s | ~10 MB |
| FMH 每个基因组（armE count） | 8 s / 基因组 | ~170–211 MB |

（db300_k16 两样本 530s vs 231s 的差异是样本本身的命中量差异，非噪声。）

## 2. 结论一：索引成本是否近同？——是，且量化如下

- **磁盘每 landmark：酶板 38.7–38.9 B，与 panel 规模无关（k16≈k8），与
  任务书的期望一致**——锚点记录就是 `(u64 hash, u32 genome, u16 contig,
  u64 pos, flags, gc)` + 12 B  packed tag，两种 landmark 源同构。
  FMH 若进入 sk2bgrow 索引格式，每 landmark 磁盘相同（~39 B）；现行
  sketch 臂的 pkl 是 23.4 B/locus，但它不存 mismatch 校验 tag、不是
  工程化二进制格式，**不构成 1.6× 的磁盘优势**。
- **构建 wall**：酶板 3.0 s/基因组（motif 扫描）对 FMH-Python 1–2 s/基因组；
  合并库 86 s / 10M anchors，峰值 RSS 3.9 GB（≈400 B/anchor 构建期，
  与 A1 的 160 B/anchor 边际量同量级口径）。FMH 阈值过滤无 motif 扫描，
  理论上在 Rust 实现里只会更快——**没有哪条路线在建库成本上输**。
- **真正的成本分野在 profile**：酶板是一个合并 DB 单次扫描（300 基因组
  230–530 s @ 2M reads）；FMH 现行实现按基因组重扫 reads（8 s/基因组，
  300 基因组 ≈ 40 CPU·min，且每个基因组要驻留一份 sketch）。工程化为
  `--mode fracminhash` 合并索引后这一分野消失；在当前实现下这是 FMH 臂
  的实际额外成本。

## 3. 结论二：半密度（scale 200）的精度代价（F5_half_density.tsv）

3 个基因组 × 4 模式（A_k16、A_k8、E_s匹配16、E_s200）× 4 深度 × 8 cell，
零失败。r / RMSE：

| 基因组 (GC) | 深度 | A_k16 | E_s匹配 | **E_s200** |
|---|---|---|---|---|
| E. coli (50.8%) | 0.5× | 0.984 / 0.365 | 0.965 / 0.454 | **0.985 / 0.371** |
| | 1× | 0.987 / 0.174 | 0.999 / 0.156 | **0.998 / 0.118** |
| | 2× | 0.999 / 0.048 | 0.999 / 0.045 | **0.998 / 0.047** |
| P. aeruginosa (65.3%) | 0.5× | 0.990 / 0.346 | 0.986 / 0.380 | **0.970 / 0.346** |
| | 1× | 0.997 / 0.140 | 0.993 / 0.146 | **0.991 / 0.112** |
| | 2× | 0.999 / 0.055 | 0.998 / 0.048 | **0.998 / 0.042** |
| S. violaceoruber (72.0%) | 0.5× | 0.995 / 0.272 | 0.965 / 0.397 | **0.961 / 0.382** |
| | 1× | 0.982 / 0.175 | 0.993 / 0.168 | **0.994 / 0.135** |
| | 2× | 1.000 / 0.027 | 0.999 / 0.045 | **0.999 / 0.040** |

- **≥1×：半密度基本无代价**（r 差 ≤0.007，RMSE 甚至更低；slope 全部
  0.89–1.03）。2×/5× 全模式 r ≥0.998。
- **0.5×：有轻微代价且随 GC 上升**：50.8% 处无损（0.985 vs 0.965 匹配档）；
  65.3% 处 −0.016（vs 酶板 −0.020）；72% 处 −0.004 vs 匹配档、−0.034 vs
  酶板。半密度把 0.5× 的边缘往失效区推近了一点，但尚未逆转结论。
- s200 相对 k16-matched 的 landmark 比中位 0.534（任务书"一半"假设成立）。

## 4. C6 预算建议

投影常量（本实验实测）：磁盘 38.8 B/anchor；GTDB species reps
136,646 × ~4.6 Mb ≈ 629 Gbp；k16 密度 9,519/Mb（300 基因组实测，
GC 24–73%）。

| 路线 | GTDB reps 磁盘投影 | 依据 |
|---|---:|---|
| 酶板 k16（A1 直系） | **~232 GB** 磁盘 | 5.99G anchors × 38.8 B |
| 酶板 k8（A6 推荐） | **~184 GB** 磁盘 | 4.75G anchors × 38.8 B |
| FMH 匹配密度（若工程化） | ~232 GB（同构） | 每 landmark 磁盘相同 |
| **FMH 半密度 s200** | **~124 GB** 磁盘 | 0.534 × 232 GB |

注意 A1 的 **752 GB 是峰值 RSS**（160 B/anchor 边际）不是磁盘；磁盘按
同口径只有 ~1/4。C6 的"内存 750 GB"若指构建期峰值，实测合并库构建
RSS ≈ 400 B/anchor（10M anchors → 3.9 GB），线性外推 GTDB 全量
~2.4 TB 构建 RSS，需分 shard 构建——这一点两种模式相同。

**建议**：

1. **匹配密度下两条路线成本等价**（磁盘同构、构建同速）；C6 选路应由
    wet-lab 与精度口径决定（F1/F2 已判：匹配密度 FMH 精度打平），
   不由成本决定。
2. **C6 的实际答案是半密度**：scale 200 把索引砍到 ~124 GB（k16 口径）
   或 ~98 GB（k8 口径半密度），而 ≥1× 精度代价 ≈ 0（3 基因组实测，
   r 损失 ≤0.007）。若产品规格保证 ≥1× 深度，半密度 sketch 是预算内
   最省的路；若必须守 0.5×，72% 高 GC 端已见 −0.03 r 的退化，需按
   目标群 GC 分布再加权（本实验只测了 3 个基因组，且 0.5× 单点差异
   的 n=8，外推到全 GTDB 需谨慎）。
3. **工程缺口**：`sk2bgrow index` 尚不接受 FMH anchor（本实验 FMH 臂的
   索引/sketch 与 profile 均走 Pilea 式 Python 管线，profile 按基因组
   重扫 reads）。要把半密度方案落成 C6 默认，需要 `--mode fracminhash`
   合并索引；这是索引格式的同构扩展，不是新数学。

## 5. 局限

- FMH 臂成本是 Python 实现的成本；Rust 工程化后构建只会更快，磁盘
  结论（同构 ~39 B/anchor）不变。
- 300 基因组为 GTDB R226 随机抽样（非 species rep 全集），密度投影用
  实测 9,519/Mb，比 A1 用的 E. coli 单点 9,422/Mb 略高，投影差 <1%。
- 半密度精度只测 3 个基因组（50.8/65.3/72.0% GC）、每 cell n=8；
  GC ≤35% 的低密度区未在半密度臂内复测（F2 已示 25–30% GC 是
  所有 landmark 模式的 0.5× 失效区，半密度只会更差）。
- profile 成本用 2M reads 单端子样本；深度线性放大的口径与 F1/F2 一致。
