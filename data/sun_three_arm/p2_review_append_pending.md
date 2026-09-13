

## 6. 重大更新（2026-09-07）：dedup 是 C 臂压扁的主因 —— 修正第 3 节 B-vs-C 的归因

来源：作业 4006836（`scripts/c_arm_diag.py`，日志 `logs/P2/cdiag.4006836.out`）。

### 结论翻转

第 3 节原来把 B-vs-C 的差归因于"真实 2bRAD 文库本身的系统性差异（建库/扩增/去重）"。
诊断把四个候选拆开（raw counts vs deduped counts × v_shape/glm × e_prior），
**raw C（不去重）vs B 远好于 deduped C vs B**：

| arm | S01 | S06 | S07 |
|---|---|---|---|
| C deduped vs B（原报告） | r 0.297 CCC 0.185 bias −0.73 | r 0.396 CCC 0.166 bias −1.13 | r 0.524 CCC 0.285 bias −0.55 |
| **C raw vs B** | r 0.449 CCC 0.440 bias −0.23 | **r 0.797 CCC 0.728 bias +0.60** | **r 0.734 CCC 0.566 bias +0.62** |
| C deduped + GLM vs B | r −0.025 | r 0.053 | r 0.412 |

口径与第 3 节相同（fraction≥0.5、dispersion<5 双边；n=8/10/16，比 deduped 臂的
18/18/20 少，因为 raw 臂过闸的种更少——这一点本身也是诊断信息）。

**awk 精确去重（95.7% 的 read 折叠到唯一序列）压扁了 PTR 动态范围，是 C 臂压扁的
主因**，不是建库偏差。旁证：σ_eff 在 deduped C 上只有 1 个种可估，raw C 有 52–116 个
——去重把计数结构整个改变了。

两个连带修正：

1. **"e_prior 首跑无改善"的归因要改**。当时写"原因明确：e_prior 只覆盖 0.3–0.7%
   锚"。这只对了一半。pooled-e 现在证明同样被门槛卡死（见下），而即便覆盖够，
   在 dedup 压扁的信号上做校正也是修错了变量——dedup 才是主因。
2. **C8 的 glm 结论在真实数据上得到独立复现**：glm 在 deduped C 上是垃圾
   （S01 r=−0.025），与 C8 收场（glm 不全局采纳）一致。glm 救不了压扁，
   压扁不是低计数斜率病理。

### pooled-e 的失败是门槛 bug，不是方法失败

pooled-e 只覆盖 9 个基因组 / 3,617 个锚：脚本对 pooled 后的表仍套用 mC≥20 门槛，
而 C 臂每样本覆盖 ~1.5×，池化 3 样本后 mC 中位 ~4.5，全被门槛杀掉。修法是
pooled 后降门槛（m≥5–8），不是放弃池化。这条留作待办，优先级低于 dedup。

### [D] bias vs coverage 未受控，暂不引用

c_arm_diag 的 [D]：Spearman(bias, log10 coverage_C) = −0.459（n=56），但最高覆盖
bin 的 bias −3.25 由 B 臂离群值驱动（PTR_B 到 1832 的基因组）。在 winsorize
（|ba_diff|≤3）重算之前，这条不能写进论文。已加入 cap dedup 敏感性作业
（4006941）一并重算。

### 论文口径影响（三档）

- **可写论文**：raw-vs-deduped 的对比本身是方法学结果——对 2bRAD 数据做
  read 级精确去重会摧毁 PTR 信号（动态范围压扁 + bias 翻转）。这直接指导
  湿实验/生信流程：sk2bGrow 路线不应在计数前做全库精确去重。
- **成立负面**：glm 不是 C 臂的解药（与 C8 相互独立的两条证据）。
- **未决**：cap dedup 敏感性（4006941，cap∈{1,2,3,5,∞}）——是否存在一个中间
  cap 同时保住 raw 的动态范围和 deduped 的离散度控制。若存在，流程建议是
  "cap-N 去重"而不是"不去重"。

### 第 5 节待办更新

1. ~~C-vs-B bias 诊断~~ → 已定位为 dedup；细化为 cap dedup 敏感性（进行中）。
2. pooled-e 重跑（修门槛 bug：pooled 后 m≥5–8）。
3. Hou（独立批次）：σ_eff 外部验证 + 规模复现。PRJCA030517。
4. F1–F7 sketch 实验（C8 已收场，解冻）。


## 7. cap dedup 敏感性结果（2026-09-08，作业 4006941，27 min）

对 deduped counts 做 C'_i = min(C_i, cap)，cap ∈ {2,3,5}，重跑 stats 与
B 对比（口径同 §6）：

| cap | S01 | S06 | S07 |
|---|---|---|---|
| 1（=现状 dedup） | r 0.297 CCC 0.185 | r 0.435 CCC 0.190 | r 0.604 CCC 0.307 |
| 2 | r 0.758 CCC 0.564 | r 0.205 CCC 0.154 | **r 0.831 CCC 0.701 bias +0.24 SD 0.386** |
| 3 | r 0.343 CCC 0.154 | r 0.266 CCC 0.099 | r 0.259 CCC 0.106 |
| 5 | r 0.320 CCC 0.195 | r 0.434 CCC 0.188 | r 0.517 CCC 0.245 |
| ∞（=raw） | r 0.449 CCC 0.440 | **r 0.797 CCC 0.728** | r 0.734 CCC 0.566 |

**结论：不存在稳健的中间 cap。** cap=2 在 S07 上取得全部 C 条件里最好
的结果（r 0.831 / CCC 0.701 / SD 0.386），但在 S06 上最差（r 0.205）；
跨样本无常性。表级 cap 变换本身也只是 capped-dedup 的粗糙代理。维持
§6 的建议：**不做计数前全库精确去重**（raw 平均最优，且三个样本方向
一致为正相关）。

**[D-winsorized]**（|ba_diff|≤3，跨样本合并）：cap1 rho=−0.223 (n=49)、
capinf rho=**+0.409** (n=34)、cap2 −0.086、cap3 −0.331、cap5 −0.221。
原始 [D] 的 −0.459 确认由 B 臂离群值驱动；winsorize 后无一致的
bias–coverage 关系，**该相关不写入论文**。
