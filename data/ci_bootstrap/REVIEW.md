# 任务 5（bootstrap CI + 2×2 交互检验）结果

日期：2026-09-11 ｜ 评审依据：甲 M3/M9 ｜ 方法：media bootstrap（16 media 重抽
10,000 次，complete-case 排除 RUN_OUT 对照），配对 bootstrap 计算 arm 间 Δr，
交互作用在 Fisher-z 尺度上检验 ｜ 脚本：`/tmp/ci_bootstrap.py`（可向论文仓迁移）
｜ 随机种子：20260911

## 判据结论

### 1. "≥1× 两 landmark 源不可区分" —— 成立，且现在可量化

F1 harness 内配对 Δr（A_k16 − E_s104，密度匹配 9,645 vs 9,422 landmarks/Mb）：

| depth | Δr point | 95% CI |
|---|---|---|
| 1× | −0.000 | [−0.163, +0.115] |
| 2× | −0.030 | [−0.112, +0.035] |
| 5× | −0.006 | [−0.042, +0.024] |
| 10× | +0.003 | [−0.022, +0.029] |

全部 CI 紧致地跨 0。"indistinguishable"从目测升级为有 CI 支撑的等价性陈述。

### 2. 0.5× 排序 —— 方向与正文一致，但 CI 很宽；跨 harness 差异不能被抽样噪声解释

- A_k16 − E_s104 @0.5×：Δ = −0.130 [−0.656, +0.257]（n=10 存活细胞）——
  sketch 在存活细胞上靠前的方向与正文一致，但 95% CI 跨 0，不能作为定论。
- A_k16 @0.5× 在 F1 harness 的 r = 0.568，95% CI [0.14, 0.89]。
  **跨 harness 差距（0.913 vs 0.568）远超 CI 宽度 → 差异是 harness 配置效应，
  不是抽样噪声。**任务 1（对称重跑）仍然必须做，CI 无法替代它。

### 3. 2×2 交互作用 —— 在论文关注的深度带内形式检验显著

interaction_z = (z_A − z_B) − (z_E − z_C_relaxed)，同一 media 重抽：

| depth | z | 95% CI | 显著 |
|---|---|---|---|
| 1× | +1.17 | [+0.30, +1.95] | 是 |
| 2× | +0.77 | [+0.09, +1.23] | 是 |
| 5× | +0.36 | [−0.32, +0.76] | 否 |
| 10× | +0.84 | [+0.20, +1.53] | 是 |

交互主张在 1–2×（论文的核心深度带）和 10× 通过检验；5× 不显著
（该深度四个格子本来就接近：A 0.981 / E 0.978 / B 0.918 / C 0.956）。
0.5× 无法计算（C_relaxed 全为常数，r 无定义）——交互的 0.5× 证据
保持"只有组合单元可用"的描述性口径。

### 4. "r 在 k=8 达峰" —— 点估计成立，统计上不可与 k=4/k=12 区分

r_low（0.5/1/2× 平均）bootstrap：

| k | r | 95% CI |
|---|---|---|
| 2 | 0.953 | [0.916, 0.977] |
| 4 | 0.960 | [0.931, 0.977] |
| 8 | 0.970 | [0.947, 0.984] |
| 12 | 0.970 | [0.943, 0.986] |
| 16 | 0.960 | [0.928, 0.979] |

Δ(8−4) = +0.010 [−0.019, +0.042]；Δ(8−16) = +0.010 [−0.019, +0.044]。
**诚实口径：4–12 酶在统计上不可区分，16 酶不带来增益，k=8 是点估计峰值
但不可与 4/12 区分。** 正文 "r peaks at k=8" 应改为 "plateau 4–12, peak at 8
within noise"。

### 5. Table 2 头条数字的 CI（arm A，参考 growth_rate）

| depth | r | 95% CI | slope | 95% CI |
|---|---|---|---|---|
| 0.5× | 0.917 | [0.832, 0.971] | 0.614 | [0.498, 0.735] |
| 1× | 0.983 | [0.960, 0.994] | 0.779 | [0.681, 0.861] |
| 2× | 0.983 | [0.957, 0.994] | 0.839 | [0.766, 0.911] |
| 5× | 0.981 | [0.946, 0.995] | 0.918 | [0.826, 1.000] |
| 10× | 0.971 | [0.941, 0.989] | 0.948 | [0.820, 1.089] |

- 低覆盖斜率压缩统计上确凿（0.5× slope CI 上限 0.735，1× 上限 0.861）。
- 5×/10× 的 slope CI 已触到 1.0——"深度处斜率趋于无压缩"的口径与 CI 兼容。

## 方法备注

- Fisher-z 渐近 p 值未采用：media 间配对且 n=16，独立性假设不成立；
  配对 bootstrap CI 是主要推断依据（fisher_z_p 列保留在 f1_ci.tsv 中，均为
  NaN 渐近值，不引用）。
- RUN_OUT 对照无实测 growth rate，complete-case 排除（n=16）。
- F1 单 arm 点估计与已提交 F1_results.tsv 逐位一致（0.9060/0.5681/0.9089/
  0.8183），确认 F1 的 r 参考系为 growth_rate。
- 0.5× 的 E_s104 仅 n=10 存活细胞，其 CI 宽主要由此而来，非脚本问题。

## 对稿件的直接修订建议

1. Table 2 加 r/slope 的 bootstrap CI 列（数据在本目录 table2_ci.tsv）。
2. §4 "indistinguishable"处加一句 CI 支撑（Δr CI 表可进补充材料）。
3. §5 "r peaks at k=8" 改为 plateau 口径（"4–12 enzymes statistically
   indistinguishable; k=8 is the point-estimate peak"）。
4. Methods §5 加 bootstrap CI 与配对 Δr 的定义段落；声明 F1–F5 为
   post-hoc、核心实验（Zheng titration / simulation / fragmentation /
   Sun 三臂）为 pre-specified（甲 M9 所需素材至此齐备）。
5. 交互检验结果（1×/2×/10× 显著，5× 不显著）进 §4 或补充材料。

---

## 写回记录（2026-09-11，同日完成）

上述修订建议 1–4 已全部写回稿件与生成器：

- Table 2/Table 7 增加 bootstrap CI 列（`figures/make_tables.py` 读
  `data/ci_bootstrap/*.tsv` 合并，表格已重新生成）；
- §2.1 头条 r 加 CI；§4 "indistinguishable" 改 CI 口径并加入交互检验结果；
- §5 峰值措辞改为 plateau 口径（4–12 酶统计不可区分，k=8 为点估计峰值）；
- Methods §5 末尾新增 "Uncertainty quantification and analysis provenance"
  段（bootstrap 定义、配对 Δr、交互检验式、pre-specified vs post-hoc 声明）。
- manuscript.md / methods.md / results_draft.md / make_tables.py 四处同源同步。

待任务 1（F1 对称重跑）返回后：若 0.5× 排序翻转，§4 的交互表述需按
"≥1× 交互 + 0.5× 融合冗余"的边界重写；CI 框架不变。

---

## 再生成与交互检验更新（2026-09-13）

Table 2 全部五臂已在 HPC C1 双端实例上用现行代码再生成（`data/repro_check/`，Mac 端旧实例备份于 `results_raw_mac_backup.tsv`）。本目录的 table2_ci.tsv 已按新 results_raw 重跑（同种子 20260911）；**interaction.tsv 是单实例（16 media）分析，其"1×/2×/10× 显著"形态已被多种子实验取代**：3 实例 × 16 media = 48 单元的配对 bootstrap 显示交互仅 2× 显著为负（z=−0.494 [−0.728, −0.278]，z_norm=−4.36，三实例同号），1×/5×/10× 不显著，0.5× 不可检验（C_relaxed 退化）。裁决报告：`data/repro_check/multiseed/INTERACTION_REPORT.md`。稿件 §4/Discussion/Methods 已按新口径重写。
