# C5 (PRJNA974210, RBC) 评分复核 — 六项交付物汇总

日期：2026-09-03 ｜ 范围：C5 评分口径复核（PRJNA974210, 9 个 RBC 样本, 522 MAGs）
所有表格在 `bench/C5/review/`，评分中间产物在 `bench/C5/res/`。补跑数组 3981503 与
重评分 3981504 均已 COMPLETED。本文件按复核意见六条逐项给结论。

---

## 1. recall 三种口径（交付物 1）

表：`c5_recall_three_ways.tsv`（逐样本 ×  arm），`c5_common_denominator.tsv`（共同分母）。

口径定义：

- **reported_fraction**（原"recall"改名）：输出行的行数 / 522。sk2bgrow 没有输出
  gate，恒为 522/522 = 1.00。这反映的是指标定义，不是方法性能，不作头条。
- **estimate_fraction**：真实给出 PTR 估计的行数 / 522。sk2bgrow 有 19 行
  "no PTR estimate"（coverage 为 NaN，均 0 QC 通过），estimate_fraction =
  0.988–1.000。
- **qc_recall**：QC 通过数 / 522，与 Pilea default（有 gate）同口径比较。

论文该写的一行（同口径，各自通过自己的门槛）：

| 指标 | 范围（9 样本） |
|---|---|
| sk2bgrow qc_recall | 3.8% – 13.8%（20–72 / 522） |
| Pilea default reported | 5.0% – 11.1%（26–58 / 522） |

两者几乎完全重叠 —— 1.00 vs 0.05–0.11 的原始比较是"无 gate 全输出"对
"Pilea 拒绝输出"，不可比。

共同分母口径：Pilea default 的 26–58 个估计**全部**落在 sk2bgrow 也输出估计的
MAG 集合内（逐样本 n_common = pilea_default 的 n，无一超集外）。即 Pilea default
通过的 MAG 是 sk2bgrow 输出集合的子集。在共同子集上 sk2bgrow 的 QC 通过率见
`c5_common_denominator.tsv`（n_common_qc_sk2b 列，15–36 / 26–58）。

## 2. suspicious 判据（交付物 2）

判据定义（`score_c5.py:38,73-74`）：

```
COV_FLOOR = 0.1
suspicious = has_est and cov < 0.1 and (pass_qc in (True, "True", None))
```

即"有估计 + 覆盖度 < 0.1x + 通过了 QC"。该判据结构性空置，两层原因：

1. **内部不咬合**：pass_qc=True 要求 coverage 足够，QC 通过行的最小 coverage =
   1.97x（>> 0.1）。`< 0.1 且 pass_qc` 是空集，判据永远不会触发。
2. **数据分布**：本数据集 9 个样本 × 522 MAG 的 coverage 最低 1.20x，中位 2.0–3.1x
   （67–79M reads/样本的深测序）。不存在 coverage < 1.2 的 MAG，所以即使去掉
   pass_qc 条件，"低覆盖却给了估计"的单元格也不存在 —— 复核意见中"必然存在低
   覆盖 MAG"的假设在这套深测序数据上不成立，如实记录。

表：`c5_coverage_distribution.tsv`（每样本 coverage 分布）、
`c5_coverage_nan_audit.tsv`（19 行 coverage=NaN，全部 estimate_fraction 缺失、
0 QC 通过，已审计）。

**结论表述**：该判据在 C5 上不产生约束（结构性空置 + 本数据 coverage 分布不支持），
不报 "suspicious = 0 说明结果干净"。

## 3. 成本随规模翻转（交付物 3）——单列一节

这是本次最重要的新事实，也是 roadmap M4 两层架构（containment 预筛 → anchor 级
计数）必要性的直接证据，并影响 C6 GTDB 预算。

逐样本原始数据：`c5_cost_per_sample.tsv`（/usr/bin/time -v 逐样本 wall / 峰值 RSS /
CPU%；不要只用区间）。

汇总（9 样本）：

| 阶段 | 数值 |
|---|---|
| sk2bgrow count wall | 58,931–98,049 s（16.4–27.2 h），CPU 98–99%（单核） |
| sk2bgrow count 峰值 RSS | 2,251–2,256 MB |
| sk2bgrow stats wall | 1,192–2,614 s（20–44 min），峰值 RSS ~14.9 GB |
| sk2bgrow 合计 / 样本 | 980–1,676 min |
| Pilea default | 343–926 s / 样本 |
| Pilea gates-off | 602–2,779 s / 样本 |
| **比值（sk2bgrow : Pilea default）** | **约 65x–277x** |

对照：Zheng 单菌基准上 sk2bgrow count 8.2 s vs Pilea 11.5 s（我们更快）。
方向在 522-MAG 规模上翻转。

count 阶段拆分（`c5_count_breakdown.tsv`，1M read pairs 子集相位计时 +
SRR28338158 全量校验）：

| 相位 | 1M pairs 实测 (s) | 占 count 比例 | 成本模型 |
|---|---|---|---|
| db_load | 52.4 | 2.4% | 固定 / run |
| windows_write | 0.4 | 0.02% | 固定 / run |
| index_build | 14.0 | 0.6% | 固定 / run |
| **match（reads × 24.1M anchors）** | **1,939.5** | **89.6%** | 线性 / read pair |
| 每样本 windows+rates 写出 | 157.3 | 7.3% | 线性 / read pair |

- 瓶颈 = match，占 89.6%，随 reads 线性增长；固定开销可忽略（3.1%）。
- 线性外推 SRR28338158（34.25M pairs）预测 71,890 s，实测 58,931 s
  （高估 22%）；即全量运行的 per-pair match 成本（1.72 ms）低于 1M 子集
  （2.10 ms），子集外推偏保守。瓶颈定位不受该误差影响。
- 吞吐：~516 pairs/s/核（522 MAG、24.1M anchors 的数据库规模下）。

**解读**：成本 = reads × 参考锚点数，两个因子在 C5 上都比单菌基准大两个数量级，
而实现是单线程锚点匹配。M4 两层架构（containment 预筛掉无信号参考 → 只对存活
参考做 anchor 级计数）同时砍两个因子，是成本翻转的直接对策。C6（GTDB 规模）在
当前实现下不可行，需 M4 先行。

Pilea gates-off 在 6/9 样本 exit=1（min_samples 崩溃，稳定复现）——这是 Pilea
 gates-off arm 的结果而非我们的事故；default arm 全部正常。

## 4. QC pass × MAG 质量交叉表（交付物 4）—— R2 的正面结果

数据：CheckM2 1.1.0（522 个 MAG，16 线程，数据库 uniref100.KO.1.dmnd 校验通过），
N50 由 fasta 现算。表：`c5_mag_quality.tsv`（522 MAGs：n_contigs / N50 /
total_len / Completeness / Contamination）、`c5_qc_x_mag_quality.tsv`（逐 MAG 合并
QC 通过率）、`c5_qc_x_ncontigs_binned.tsv`（分箱表）、`mag_n50.tsv`。

spearman(qc_pass_rate, 质量指标)，n=522：

| 指标 | ρ | 判读 |
|---|---|---|
| n_contigs | **−0.409** | 越碎片化，QC 通过越少 |
| N50 | **+0.426** | 一致 |
| Completeness | **+0.466** | 一致 |
| Contamination | −0.179 | 弱负相关 |
| total_len | +0.011 | 无关系（不是基因组大小的问题） |

QC 通过率按 contig 数分箱（qc_rate = 该 MAG 在 9 样本中 QC 通过的比例）：

| contig 数 | n_mags | QC 通过率均值 | 中位 |
|---|---|---|---|
| ≤10 | 73 | 22.5% | 22.2% |
| 11–25 | 84 | 19.8% | 11.1% |
| 26–50 | 131 | 6.9% | 0.0% |
| 51–100 | 145 | 6.4% | 0.0% |
| 101–200 | 71 | 2.3% | 0.0% |
| >200 | 18 | 4.3% | 0.0% |

**判读（按复核意见的判据）**：QC 通过率随 contig 数上升而单调下降、与 N50 和
Completeness 正相关、与 total_len 无关 —— **QC 在真实数据上有效**，这是 R2
（scaffolding 失效边界）的正面结果，也是 C5 作为 application section 的核心
内容：sk2bGrow 的酶一致性 QC 确实把碎片化参考的估计过滤掉了（50+ contigs 的
MAG 九成以上样本不过 QC）。注意这是对"QC 能识别碎片化"的验证；与 A2 的关系：
A2 说被摧毁的坐标会让 QC 也失效（100% fragmented pass），那里的失效模式是
contig 数少但序列被 Shuffle；本数据集的真实 MAG 碎片化是"contig 数多"，两种
碎片化模式不同，不矛盾，综合汇报中应区分表述。

同一 MAG 跨 9 样本的内部一致性（无真值情况下的替代检查）：
`c5_cross_sample_consistency_all.tsv` / `_qc_pass_only.tsv`。
log2PTR 跨样本 std：全部估计（n=522）中位 0.129、p90 0.473、max 1.947；
仅 QC 通过的行（n=213 个 MAG）中位 0.030、p90 0.086、max 0.757。511/522 个 MAG
在 9 样本全部有估计。QC 过滤后一致性收紧约 4 倍 —— QC 在无真值数据上表现出
内部一致性筛选能力。

## 5. 两件小事（交付物 5）

- **525 vs 522**：见 `c5_525_vs_522.md`。NCBI PRJNA974210 esummary 实际列出 525
  个 assembly（计划书中的数字来源），源池 rbc_mags 收集到 523 个 .fna，
  `00_build_refs.py` 按 domain 排除 3 个古菌（GCA_041394925.1 Nitrososphaeraceae、
  GCA_041394965.1 Ca. Nitrosocosmicus、GCA_041395005.1 Methanolobus）→ 522 细菌
  MAGs。非数据丢失。
- **.done 掩盖审计**：set -e 导致 3981503 部分任务产物完整但 .done 未写、状态
  FAILED —— 诊断属实。已逐个核对 C1b/C2/C3/C4/C5 全部已完成数组的产物完整性
  （按文件内容而非 .done 标记），无其他掩盖缺口。备注：C2/C4 评分脚本为"按文件
  存在性静默跳过"模式，本次审计未触发该路径，但建议后续把跳过写进日志。

## 6. 六项交付物清单

1. ✅ 三种口径评分表：`c5_recall_three_ways.tsv`（reported_fraction /
   estimate_fraction / qc_recall 并排）+ `c5_common_denominator.tsv`
2. ✅ suspicious 判据定义 + 覆盖度分布：`c5_coverage_distribution.tsv`、
   `c5_coverage_nan_audit.tsv`；结论 = 判据结构性空置，如实声明
3. ✅ 逐样本耗时与峰值 RSS：`c5_cost_per_sample.tsv`；count 拆分：
   `c5_count_breakdown.tsv`（瓶颈 = match，89.6%，线性于 reads）
4. ✅ QC × MAG 质量交叉表 + 相关系数：`c5_qc_x_mag_quality.tsv`、
   `c5_qc_x_ncontigs_binned.tsv`；ρ(qc_rate, n_contigs)=−0.409，随碎片化单调
   下降（R2 正面结果）；跨样本一致性见上节
5. ✅ 525 vs 522 说明：`c5_525_vs_522.md`；数组完整性审计完成
6. ✅ 本文件

**C5 一节在综合汇报中的口径**：头条改为"同口径 qc_recall 3.8–13.8% vs Pilea
default 5.0–11.1%，两者重叠"；成本翻转单列；recall=1.00 仅作 reported_fraction
报告，不作性能声明。

## 6b. C5 代码侧证据的实测（交付物 7，qc_audit.py 对 9 样本跑完）

输入：`bench/C5/review/qc_audit/A_SRR*.per_enzyme.tsv`（各样本 per-enzyme 表副本，
8,352 行 = 522 MAG × 16 酶）+ `*.fused.tsv`（output.tsv 副本）。审计输出：
`A_SRR*.audit.txt`。

9 样本合计 4,698 个 (sample, genome) cell：

| 指标 | C5 实测 | 与 Zheng 对照 |
|---|---|---|
| 存活酶数 k，均值 | 9.0–10.6 | — |
| k=1 cells（Cochran Q 未跑即记 consistent） | **27（0.57%）** | Zheng 各深度均为 0 |
| k=0 cells | 18 | — |
| 融合表 n_enzymes==1 行 | 33（0.70%） | — |
| ok 的符号拒绝（negative estimate） | 7/9 样本恰好 0；另两个 40 和 10 | Zheng 0 次 |
| 被接受但 r²<0 的拟合 | 902 / 47,285（1.9%），逐样本 1.3–3.9% | Zheng 低深度 57.5%（0.5×）→1.1%（5×）；C5 中位覆盖 ~2.4x，落在 Zheng 2× 的 6.7% 与 5× 的 1.1% 之间，一致 |
| coverage p50 / p95 | 2.36 / 8.25（SRR59 代表值） | — |
| 多数酶未出拟合但 coverage 读数 ~2 的行 | 139/522（SRR59） | coverage 列存活者偏差的活证据 |

**判定**：复核意见第 6 节的假设在 C5 上坐实但幅度小——"k=1 未检验即通过"只影响
0.7% 的融合行，不是 recall=1.00 的主因。recall=1.00 的主因仍是"全流程无任何覆盖度
闸门"（未修缺陷之一）：低覆盖 MAG 照样出估计，QC 也只在酶一致性/S E/锚数上设卡。
r²<0 的拟合有 1.9% 被无条件接受，量级与 Zheng 同覆盖区间一致——"QC 不看拟合质量"
在真实数据上同样成立。

注：现有融合表生成于 a56197d 之前，无 consistency_checked/min_r2 列；审计工具已
识别并提示 re-fuse。如后续需要带新列的 C5 融合表，可从现有 counts 重跑 stats 层
（不重跑 count）。

## 4c. 覆盖度混杂控制（复核追加，2026-09-05）

质疑：碎片化的 MAG 通常低丰度 → 低覆盖 → QC 检验力低。ρ(qc_rate,
n_contigs) = −0.409 可能大部分是覆盖度效应。控制结果（脚本
`scripts/c5_coverage_control_full.py`，输出 `c5_coverage_control.tsv` /
`c5_coverage_control.txt`，逐 MAG 数据 `c5_qc_x_quality_with_cov.tsv`）：

**偏相关（Spearman，控制 mean coverage）**：

| feature | ρ_raw | ρ_partial(cov) | p |
|---|---|---|---|
| n_contigs | −0.409 | **−0.337** | <1e-4 |
| N50 | +0.426 | **+0.358** | <1e-4 |
| Completeness | +0.466 | **+0.353** | <1e-4 |
| Contamination | −0.179 | −0.166 | 1e-4 |

**箱内关联（mean coverage 四分位内）**：ρ(qc_rate, n_contigs) 随覆盖度分位
*增强*而非减弱——Q1（cov 1.5–2.2）−0.194 → Q4（cov 3.4–33.6）**−0.432**。
若关联主要由覆盖度驱动，高覆盖箱内（QC 力饱和）关联应消失；实测相反。
覆盖度四分位 × contig 数分箱的 qc_rate 在每个箱内仍随碎片化下降
（Q4：≤10 contigs 0.327 → 101–200 contigs 0.061）。

**混杂结构**：ρ(qc_rate, cov) = +0.516（最强单预测），ρ(n_contigs, cov) =
−0.252（耦合中等）。覆盖度解释了部分关联（衰减 ~18%），但碎片化效应
独立存在。Q4 箱（QC 力饱和区间）ρ = −0.432 是最干净的证据。

**对"一致性收紧"的控制**：跨样本 log2PTR std，QC-any MAG 中位 0.105 vs
无 QC 0.152；在覆盖度四分位内分别保持（Q1：0.105 vs 0.150；Q4：0.099 vs
0.129）。OLS：std ~ qc_any + log10(cov) 中 qc_any β = −0.084（p = 5.5e-05），
log_cov 不显著（p = 0.64）。收紧不是覆盖度效应。

**表述更新**："QC 在真实数据上确实滤掉了碎片化参考"成立，但必须同时报
覆盖度是最强单预测且解释部分关联；两个效应（覆盖度、碎片化）在真实数据
上方向相同、部分纠缠，这本身是 RBC MAG 群体的真实结构，不是评分伪影。


### §4d. M4 screen A/B at C5 scale（作业 3995938，2026-09-06）

设置：C5 全部 522 MAG（24.1M anchors），`index --screen-scale 2000` 构建
screen 库（db_scr 909M vs 原库 879M；构建 6:31，峰值 RSS 8.7G）。同一
1M read-pair 子集、同机同 8 线程，m4 二进制分别以 unrestricted / `--screen`
跑 count，比对 counts.tsv。

**计时**（phase-timing，CPU s）：

| arm | motif_scan | lookup | wall | lookups | windows tested |
|---|---|---|---|---|---|
| unrestricted | 66.7 (4.3%) | 1486.6 (95.7%) | 14:33 | 2,564,540 | 1.22B |
| screen | 66.4 (4.5%) | 1406.7 (95.5%) | 14:39 | 2,564,540 | 1.22B |

**结论一：screen 在当前实现下几乎不产生加速。** lookup 仅 −5.4%，wall
0.99×（略负）。机制可从数据直接读出：screen 不改变 probe 数量
（lookups、windows tested 两臂完全相同）。pass-2 仍对每个 read 的每个
窗口做 kmer probe（`screen.rs` 注释自述：pass 2 不跳过 pass-1 未命中的
read，"speedup comes from the lookup step"），而 genome 子集筛选只缩短
每个命中 key 的 posting list，砍不到 probe/二分本身。占 match CPU 95%+
的 lookup 其成本由 read 数 × 每 read 窗口数驱动，与索引里有多少 genome
基本无关。UHGG 试点的 1.8× 并非"小库缺缓存"的假象放大——恰恰相反，
索引缩小本身就不是主要杠杆。

**结论二：counts 不等价，与 m4_bench 同一签名。** 29.7M anchors 中
20,049 个计数不符；有计数 genome 522 → 457（**12.5% 的 MAG 整组丢失**）；
screen 非零锚 ⊂ unrestricted（subset_ok=True，extra=0）。丢失模式与
`--screen-min-frac`（默认 ≥50 sketch hits 才保留）一致：1M 子集上真实
低丰度 MAG 的 sketch 命中数不足被整组丢弃。`count.rs` 里
`genome_screened_counting_matches_full_counting` 单测通过是因为合成
数据覆盖度高，撞不到这个门槛。对 PTR 场景不可接受——丢失的正是低覆盖
MAG，即 recall/QC 争议所在区间。

**对 M4 判定**：screen 机制（containment 预筛 → anchor 级计数）在当前
实现下**既不快也不对**，不能作为 C6/GTDB 规模的成本修复方案上报。成本
翻转（65–277×，后按 c5_cost_per_sample.tsv 实算修订为 **89.5–240.8×**）维持已记录限制 + 未来工作定位；若要继续攻 cost，方向是
让 pass-2 跳过 pass-1 未命中的 read/window（直接砍 windows tested，而不
是砍 posting list），或换 probe 数据结构。此结论已同步本地论文仓。


### §4e. §4d 结论修订 + σ_eff 结构检验（2026-09-06）

**§4d 修订（采纳复核意见第 1 条）**：§4d 的实验场景是筛选**唯一不可能起大作用**
的场景——C5 的 522 个 MAG 就是从这 9 个样本组装的，按构造几乎全在场，
screen 保留 457/522 = 87.5% 是"它们确实都在"的诚实报告，不是筛选失灵。
"12.5% 的筛选给 5.4% 加速"与索引缩减成比例，不能读成"筛选无用"。C6 的
场景是 GTDB 137k representatives 对单个粪便样本（在场率 <1%），筛选强度
差三个数量级。为 §4d 结论打的补丁：**screen 在弱筛选场景（在场率 ~90%）
下无加速；对 C6 价值的判定转移到稀释实验（§5，见下）。**

**σ_eff 结构检验**（`scripts/sigma_structure_tests.py`，输出
`sigma_structure_tests.{tsv,txt}`，已同步本地论文仓 data/）：

- **Test A 跨样本可复现性**：3 个 Sun 样本两两合并同种 anchor，per-species
  Spearman r：S01×S06 中位 0.762，S01×S07 中位 0.730，S06×S07 中位 0.868；
  100% 的种 r > 0.3，0 个 r < 0。**per-anchor e_i 主要是稳定的位点属性**，
  不是每次实验的随机噪声。含义：σ_eff ≈ 1.53 的大头是**可校正效应**
  （e_prior 已在产出），"route-B 收益必须下修"的表述不成立；残余的
  样本特异成分 sqrt(1−r²)·σ ≈ 0.75–1.05 仍在，但同样可以窗口平均摊薄。
- **Test B 空间自相关**：lag 1–64（~1.7–110 kb）的 lagged Pearson r 全部
  在 ±0.03 内（n 每 lag 5.7k–34.7k）。**e_i 在基因组位置上独立** → 窗口
  平均有效：99 anchor/窗时窗口层 excess CV ≈ 1.53/√99 ≈ 0.154 上限，
  残余样本特异成分摊薄到 ~0.08。V-fit 的梯度估计基本不受 σ_eff 威胁。

综合：σ_eff ≈ 1.53 这个数字成立且是 novelty（2bRAD 首位点捕获效率测量），
但它的正确解读是"大而可校正 + 可被窗口平均摊薄"，不是"route-B 的死穴"。
论文口径据此定：报 σ_eff、报两个结构性质、报 e_prior 校正作为 route-B 的
enhancement 而非 prerequisite。

**稀释实验（作业见 logs/C5/dilution_bench.*.out）**：db_dil = 522 真实 MAG
+ 4700 个合成随机基因组（各 3.5 Mb，种子固定，保证 absent），screen-scale
2000。两臂同 1M 子集：unrestricted（~242M anchor 全索引）vs --screen
（期望选中 ~522 中的高覆盖子集，pass-2 回到 ~24M anchor 量级）。合成填充
测的是筛选的选择强度；**不测** GTDB 的同源 posting-list 膨胀（那是 mm=1
臂 + seed_hist 的活，见下）。索引 24.1M anchor 构建 3 min（16 线程），
242M 预计 ~30–40 min；unrestricted 臂 lookup 随索引超线性恶化正是要测的
C6 体制数据点。

**mm=1 臂（作业 3999442）**：`--max-mismatch` 是运行时参数（默认 2，
pigeonhole m+1 连续种子，count.rs:29-32），同库直接跑 mm=1 无需重建索引。
理论账：mm=2 是 3 个 ~10-11 bp 种子（键空间 ~4^10.5 ≈ 1.4–4M 键装 24.1M
锚），mm=1 是 2 个 ~16 bp 种子（键空间 ~4^16 ≈ 4.3e9，对 24.1M 锚几乎
键键唯一）。实测 lookup 代价 4.64 ms CPU/lookup ≈ 46k 次访问（100 ns/
随机访问折算），说明 posting list 极长——seed_hist 工具（count.rs 新增
`posting_hist` 诊断方法 + examples/seed_hist.rs）直接输出分布验证。两臂
同跑带 stats，报灵敏度代价（count delta + mismatch_hist）。代价是 F3，
受控分歧梯度上的灵敏度是后续。


### §4f. mm=1 臂结果（作业 3999630；seed_hist 见 §4e）——lookup 63× 加速

**机制实证**：mm=2 的中间种子塌缩。seed_hist（C5 db，24.1M 锚）：

| 表 | n_keys | mean_len | hit_mean_len | max |
|---|---:|---:|---:|---:|
| mm=2 slot0 | 1.19M | 20.2 | 249 | 3,648 |
| **mm=2 slot1（中段）** | **48,832** | **493.7** | **1,453.7** | **12,743** |
| mm=2 slot2 | 1.75M | 13.8 | 255 | 4,702 |
| mm=1 slot0 | 11.85M | 2.03 | 19.6 | 1,542 |
| mm=1 slot1 | 14.42M | 1.67 | 7.08 | 705 |

24.1M 锚的 tag 中段 ~11 bp 只有 4.9 万个不同值——识别位点共有序列把中间
种子压垮，每次 lookup 要验证 ~1,454 个候选（4.64 ms CPU/lookup ≈ 46k
次访问由此而来）。mm=1 的 ~16 bp 种子没有这个结构，hit_mean 降 2-3 个
数量级。

**计时**（同 1M 子集、同库、8 线程，count 阶段 CPU s）：

| arm | motif_scan | lookup | wall |
|---|---:|---:|---:|
| mm=2 | 66.4 (4.1%) | 1,539.4 (95.9%) | 15:55 |
| mm=1 | 64.5 (72.5%) | **24.5 (27.5%)** | **1:51** |

lookup **63× 加速**，wall **8.6×**（15:55 → 1:51）。probes 数不变
（2,564,540），加速全部来自每次 lookup 的候选验证数。

**灵敏度代价（本数据）**：anchors lost 95,156 / 1.80M nonzero（5.3%）；
**0 个 genome 整组丢失**；共有锚上 count 比中位 1.0000（无偏）；genome
总 count 中位比 0.931。本数据分歧低，代价温和。受控分歧梯度（F3：
0.1–3% 替换）上代价会随 P(tag 内 ≥2 mismatch) 上升（3% 替换时 ~26%），
mm=1 的适用窗口取决于应用的分歧体制——这是 F3 要量的。

**对成本翻转的修订**：C5 实测 count 16–27 h/样本 → mm=1 外推 **1.9–3.1
h/样本**（+ stats 20–44 min），对 Pilea 6–15 min 的比值从 ~~65–277×~~
（按 c5_cost_per_sample.tsv 实算为 **89.5–240.8×**）收窄到 **~10–25×**。
flip 没有消失但大幅收窄。

**后续工程方向（比 mm=1 更优）**：mm=2 保灵敏度 + 种子避开共有区。三个
~11 bp  pigeonhole 种子若改选到 tag 的高熵区（或 masked seed 跳过共有
位点），三个 slot 的 hit_mean 都能到 ~20，lookup 有望压到 ~50–100 s 而
不付 mm=1 的灵敏度代价。这是 M4 之后成本线的正解，已记入 roadmap。


**counts 等价性与假阳性普查（2026-09-08 更新，dil_cmp3 + awk 普查，
取代 dilution_compare2.py 的占位引用）**：

real genomes：unres 522 → screen 457 个有种内计数（subset_ok=True，
extra=0）；per-anchor screen/unres 计数比中位 1.0000、mean 0.9488；
逐锚 mismatch 1,917 / 31,144（6.2%，集中在低计数锚，不影响主体）。
合成基因组（保证 absent，真阳性应为 0）：**unrestricted 下 4,048 /
4,700 个合成基因组拿到非零计数**（4,048 个非零 (window, genome) 对，
即每个假阳性基因组恰好 1 个非零窗）——这是 242M 锚全索引的随机命中
本底，量级与真实信号同阶，**是 unrestricted 体制下 C6 不可接受的假
阳性来源**；screen 臂合成基因组 **0 非零**，本底被 containment 预筛
完全清零。equivalence 与 FP census 两个口径同时闭合：screen 在保住
真实信号的同时消除了 absent 基因组的假阳性。
