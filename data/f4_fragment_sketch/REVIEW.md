# F4 REVIEW — 碎片化 × sketch 模式：scaffold 对 landmark 来源是否无感

日期：2026-09-08。任务：HPC_TASKS_SKETCH_MODE.md §F4。产物：
`F4_collapse.tsv`（r/RMSE/slope 按 基因组×臂×条件×深度）、
`F4_scaffold.tsv`（放置精度 + PTR 恢复）、`F4_all_cells.tsv`（288 cell 逐条）、
`F4_permutation.tsv`（frag 机制置换检验）。

作业链（sbatch，全部 COMPLETED ExitCode 0:0）：
prep 4010662→4010668（fragment.py 碎片化 + 合成近缘）/
scaffold+index 4010680→4010690→4010701 / counts 4010706→4010712 /
stats 4010723 / aggregate 4010760 / permutation 4010774。

## 0. 口径

1. **基因组**（F2 genomes，k16 匹配 scale）：g08 E. coli 50.8% GC s106、
   g06 B. subtilis 43.5% s120、g13 P. aeruginosa 65.3% s103。
   F2 基因组是 3–4 contig，先按文件序拼接成单条伪染色体再碎片化
   （`f4_concat.py`；索引的全局坐标系统就是同一拼接，layout 真值与
   参考坐标一致）。P. aeruginosa 在 F2 集内无姐妹种，scafRel 参考用
   合成 0.1% 替换突变体（`f4_mutate.py`，7,534 SNP，seed 4，明确标注
   合成）——正是 Syn2b 的 0.1% 分歧框架。
2. **碎片化**：A2 同法（`src/benches/fragmentation/fragment.py`，100 个
   lognormal contig、乱序、独立 50% 翻转、seed 0），布局真值记 layout.tsv。
3. **reads**：复用 F2 的模拟 FASTQ（只读；planted log2PTR ∈ {0.5, 1.0, 2.0}
   × 2 seed × {1, 5}×）。每 (基因组, 条件, 深度) n=6。
4. **两臂**：酶臂 = `sk2bgrow profile --no-stats` + `--windows`（A2/C1 parity）；
   sketch 臂 = `armE_counts.py`（Pilea 解释器，k31 FMH，匹配 scale，
   cache 复用/新建）不带 `--windows`（armE/F2 parity）。
   complete 条件的酶索引直接复用 F2 db（只读）。
5. **sketch 臂 scaffold**：发布的 `sk2bgrow scaffold` 只消化酶标签，FMH 版
   放置算法按 `scaffold.rs` 1:1 移植到 Pilea FMH key 上
   （`scripts/f4_scaffold_sketch.py`：共享 tag 映射丢弃多拷贝、Kendall
   投票定向、中位数起点、min_tags=3 / min_concordance=0.8），输出同一
   JSON schema，下游 rescaffold.py 评分与重建 FASTA 完全不变。移植的
   正确性由 self 臂近乎完美的放置（98–100%，Spearman 1.0，中位误差 0 bp）佐证。

## 1. 问题一：V-fit 在 contig 上是否同样塌陷？——塌陷，且与 landmark 类型无关 ✓

**坐标信号确实消失**：frag 条件下每个酶的单酶 V-fit r² 全部为负
（complete 时 r²≈0.9）；置换检验（打乱窗口位置、用同一 `fit_v_shape`
重拟合 200 次）显示 frag 上的估计值与置换零分布同阶（g08 5× 中位：
酶臂 planted 0.5/1.0/2.0 → 实测 0.20/0.37/0.73，零分布 0.11/0.16/0.29；
sketch 臂 0.16/0.34/0.73 vs 0.27/0.26/0.24）——即 frag 上的"信号"主体是
ori 网格搜索的 winner's-curse 假 V，其幅度随窗口速率散布（∝ planted
PTR）放大，并非恢复的坐标信息。

**但 r 对 planted truth 保持 ≈1.0 是假象**：planted 值域宽（0.5–2.0）+
假 V 与 planted 成正比 → 秩相关存活。诚实的损伤指标是 **slope 与 RMSE**
（F4_collapse.tsv，5×，三基因组两臂一致）：

| 条件 | slope（complete≈1.0） | RMSE log2（complete≈0.02） |
|---|---|---|
| frag 酶臂 | 0.26–0.37 | 0.84–0.87 |
| frag sketch 臂 | 0.34–0.38 | 0.84–0.85 |

即 frag 把 slope 压到 ~1/3、RMSE 放大 ~40 倍，**两个臂数字相同**——
机制是坐标丢失，与 landmark 类型无关，符合预期。与 C2（真实 Zheng 数据）
对照：C2 中 frag100 的 r 从 0.95 掉到 ~0.0–0.75、RMSE 0.07→0.92——同一
塌陷在 λ 动态范围窄的真实数据 regime 下表现为 r 崩溃；F4 的 planted 设计
只是让 winner's-curse 垃圾看起来单调。**论文表述应沿用 A2 的"塌陷"结论，
但机制描述要精确：坐标丢失杀死的是坐标拟合（r²<0），残留估计值是搜索伪影。**

## 2. 问题二：scaffold 能否用 FracMinHash landmarks 放置 contig？——能 ✓

**self 臂**（frag 对同基因组 scaffold，F4_scaffold.tsv）：

| 基因组 | 模式 | 放置 | 定向正确 | 顺序 Spearman | 5× RMSE |
|---|---|---|---|---|---|
| E. coli | 酶 | 98/100 (99.0% bp) | 100% | 1.0000 | 0.027 |
| E. coli | sketch | 98/100 (99.0% bp) | 100% | 1.0000 | 0.033 |
| B. subtilis | 酶 | 100/100 | 100% | 1.0000 | 0.017 |
| B. subtilis | sketch | 100/100 | 100% | 1.0000 | 0.017 |
| P. aeruginosa | 酶 | 100/100 | 100% | 1.0000 | 0.030 |
| P. aeruginosa | sketch | 100/100 | 100% | 1.0000 | 0.030 |

self scaffold 后 PTR 恢复到 complete 水平（slope ≈1.0、RMSE ≈0.02，
对照 frag 的 ~0.85）。**sketch landmark 的放置与酶标签逐项打平。**

**rel 臂**（对近缘基因组 scaffold）：
- P. aeruginosa vs 0.1% 合成近缘：两臂均 100/100 放置、100% 定向、
  Spearman 1.0，PTR 完全恢复（5× RMSE 0.013/0.030 vs complete 0.020/0.030）。
  Syn2b 的"0.1% 分歧下 FMH landmark 保留率更高（94.6% vs 89.5%）"在本
  设置下未转化为可分辨的放置率差异（两端都是 100%）。
- E. coli vs E. sp005843885、B. subtilis vs B. atrophaeus（真实远缘、
  大规模重排）：放置率下降（E. coli 81/77 of 100；Bacillus 44/45 of 100），
  定向"正确率"≈0–6% 是相对于 draft 真值的度量——放置其实跟随**参考的**
  共线性框架，重排断点处的 contig 因 concordance 不足被弃置（符合
  scaffold.rs 设计）。PTR 恢复部分有效：5× RMSE 酶臂 0.149/0.644、
  sketch 臂 0.106/0.223（E. coli/Bacillus），远好于 frag（~0.85），
  不及 self/complete。**sketch 臂在 rel 情形 RMSE 略优于酶臂**
  （0.106 vs 0.149；0.223 vs 0.644），方向与任务文档"possibly better"
  一致，但 n=6 cell，只宜作方向性引用。

## 3. 结论

1. **scaffold 对 landmark 来源无感，在 sketch 模式下同样工作**（self 臂
   与酶臂逐项打平，rel 臂不差于酶臂）。论文 Discussion 里
   "scaffolding 是保留酶板的理由"在 route A **不再成立**——
   scaffolding 需要的只是 draft 与参考共享 landmark，FMH k-mer 与
   酶标签同等地满足。
2. V-fit 塌陷在 sketch 臂与酶臂相同（坐标机制，landmark 无关），
   且"塌陷"的确切含义应写为：坐标拟合死亡（单酶 r²<0），残留正估计值
   是 ori 搜索伪影（置换零分布证实），slope 衰减 ~1/3、RMSE ×40。
3. 边界：rel 臂只在参考与 draft 共线性保持良好时完全恢复 PTR
   （0.1% 分歧：是；属间/重排远缘：部分恢复）——这与 landmark 类型无关，
   是 A2 已报的 scafRel 行为的复现。

## 4. 复现

`bench/F4/scripts/`：f4_01_prep.sh（拼接+碎片化+合成近缘）→
f4_02_scaffold_index.sh（酶臂 `$BIN scaffold` + sketch 臂
`f4_scaffold_sketch.py` + rescaffold.py 评分/重建 + F4 索引）→
f4_03_counts.sh（288 counts）→ f4_04_stats.sh（288 stats）→
f4_05_aggregate.sh / f4_06_permutation.sh。全部幂等；F1/F2 产物只读未改。
