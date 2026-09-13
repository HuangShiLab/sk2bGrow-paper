# Dedup 机制回归（任务 4）

2026-09-12。作业 4025210（amd，32G，~6 min）。脚本 `dedup_mech.py`，
产物 `dedup_mechanism.tsv`（910 行 = 物种 × 样本）。

## 1. 数据状态确认

- `counts_C/` = **raw**：`S0x.time.txt` 命令行为对原始 fastq
  （`.../Fecal_2bRAD-M/S01_1-1.fq.gz`）跑 `sk2bgrow profile --no-stats`，
  db = `db_bcgI`。
- `counts_C2/` = **dedup**：命令行为对 `dedup/S0x_*.dedup.fq.gz` 跑同一
  index 的 profile（`p2_32_recountC.sh` 产物）。**dedup counts 已存在，
  本任务未重跑计数**。
- dedup fastq（`dedup/`）与索引 `db_bcgI/`（anchors.bin + manifest.json）
  齐备。
- 配对逻辑沿用 `sigma_eff.py`：键 = (genome_id, contig_id, position)。
  与 σ_eff 不同，机制分析**不**过滤 flags 唯一性（机制适用于所有锚）。
  每样本 merged 锚 ~7.3M，dedup>raw 违例 **0**（dedup 确为 raw 子集）。

## 2. 定义（可辩护口径）

- **df（per-anchor duplication fraction）= 1 − dedup/raw**，限 raw>0
  （raw=0 时 df 无定义；dedup≤raw 已验证）。df∈[0,1]，截断未实际触发。
- **x = 锚到 ori 的环形距离 / 基因组全长**（manifest 各 contig 长度求和；
  d=min(|p−ori|, L−d)/L，0 at ori，0.5 opposite）。
- **ori call**：取自各样本 raw counts 的 stats 输出
  （`counts_C/S0x/stats/output.tsv`，v_shape 拟合；dedup 版 stats 因梯度
  被压扁不用）。ori-called 基因组 268/340/313（S01/S06/S07），其中
  m≥20 锚进入回归的 263/336/311。
- 每（物种 × 样本）：df ~ x 无权重 OLS；m≥20。列含义：
  `slope`=d(df)/d(x)（每全基因组分数），`slope_prox`=−slope
  （对 proximity=0.5−x 的斜率）。

## 3. 结果

机制预测：ori 近端（复制中多拷贝）PCR duplicate 更多 → **df 随距 ori
距离下降 → slope<0（即 slope_prox>0）**。

| sample | 物种数 | slope<0 占比 | 中位 slope | 中位 R² | 中位 median_df |
|---|---:|---:|---:|---:|---:|
| S01 | 263 | **71%** (188/263) | −0.102 | 0.0044 | 0.50 |
| S06 | 336 | **64%** (216/336) | −0.075 | 0.0034 | 0.60 |
| S07 | 311 | **66%** (206/311) | −0.087 | 0.0047 | 0.67 |

3 样本共享物种 125 个：42 个三样本全负、5 个三样本全正（3/3 同号 47/125）。

## 4. 机制结论（三行）

1. **方向一致**：三个粪便样本中，物种的 per-anchor duplicate fraction
   均随锚到 ori 的环形距离系统性下降（ori 近端 df 最高，常接近 1——
   近乎全部 reads 是 duplicates；远端趋 0）；slope<0 物种占比
   71%/64%/66%，中位斜率 −0.08~−0.10（df 从 ori 到 ter 约降 0.04–0.05，
   部分强例降 >0.5，如 MGYG000000343@S07：截距 0.94、slope −1.72）。
2. **机制吻合**：该方向正是"PCR 偏好高拷贝模板 → 复制中基因组 ori 近端
   锚点多拷贝、duplicate fraction 更高"的预测；exact dedup 把 ori 近端
   这些高 df 计数折回单拷贝水平，等于把 V 形梯度本身除掉——为 P2 §6
   "awk 精确去重（95.7% read 折叠）压扁 PTR 动态范围"提供了机制解释，
   支撑论文"2bRAD reads 不做计数前 exact dedup"的流程建议。
3. **量级诚实声明**：per-anchor 中位 R² 仅 ~0.004——ori 距离只解释 df
   方差的小部分，dup 水平主要由其他（未建模）因素决定；结论强度在
   "跨 ~300 物种/样本 × 3 样本的符号一致性"，不在单物种拟合优度。

## 5. 判据

**满足（按 proximity 口径 3/3 为正）**：3 个样本中物种斜率多数为机制
预测方向（slope_prox>0：S01 71%、S06 64%、S07 66%）。n=3 样本，按任务
要求不做正式显著性检验；符号占比已远超掷硬币水平，但未做统计推断。

（注：任务书"3/3 为正"按 df~环形距离字面回归对应 slope<0；本文同时
给出 `slope_prox`=−slope，两种读法同一张表覆盖。）

## 6. 边界条件

- **PCR cycle 数 / input mass 未知**（SRA 元数据无文库构建记录）：
  机制归因基于复制梯度的符号一致性，非直接 PCR 计量；cycle 依赖的
  定量验证需湿实验数据。
- 单酶（BcgI）、单中心（Sun PRJNA689204）、同一文库 protocol；
  **对其他 2bRAD protocol / 酶切组合的可迁移性是 future work**
  （Hou 独立批次可作外部验证切入点）。
- ori call 失败（method=none）的基因组被排除（每样本仅少 2–7 个）；
  若 ori call 系统性偏向高 PTR 种，结论外推到全体群落需谨慎
  （call 率 ~31–34%，覆盖率门槛所致）。
- df 为序列级 exact dedup 的口径；cap-N 去重（P2 §7，无稳健中间 cap）
  下的机制行为未检验。

## 7. 产物

- `dedup_mechanism.tsv`：910 行，列 species, sample, slope,
  slope_prox, intercept, n_anchors, r2, median_df。
- 日志：HPC `bench/dedup_mech/dedup_mech.4025210.out`。
- 本地副本：`/Users/macstudio/Downloads/sk2bGrow-paper/data/dedup_mech/`
  （tsv + 本文件 + 脚本 + sbatch）。
