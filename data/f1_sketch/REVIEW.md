# F1 REVIEW — Density-matched FracMinHash vs enzyme panel, Zheng grid

日期：2026-09-08。数据：`F1_results.tsv`（10 arms × 5 depths 聚合）、
`F1_all_cells.tsv`（逐 cell）、`F1_mechanism{,_cells}.tsv`（机制列，
来自 windows.rates.tsv）、`out/<arm>_<medium>_<cov>x/`（stats 全产物）。
作业链：index 3989009 / counts 3989010 / sketch 3989011 / stats 3989012+
3989239 / agg 3989013+。

## 1. 结果（匹配对：A_k2↔E_s268, A_k4↔E_s184, A_k8↔E_s123, A_k16↔E_s104）

| 深度 | panel r（n） | sketch r（n），匹配密度 | 判定 |
|---|---|---|---|
| 0.5× | 0.50–0.57（15–16） | **0.82–0.84（10–11）** | sketch 条件精度更高，但丢 6/17 cell |
| 1× | 0.89–0.91（15–16） | 0.85–0.92（16–17） | 打平 |
| 2× | 0.93–0.98（15–16） | 0.92–0.97（15–17） | 打平（k2/k8 panel 略前） |
| 5× | 0.98（16） | 0.98（16–17） | 打平 |
| 10× | 0.97–0.99（16） | 0.97–0.99（15–16） | 打平 |

landmarks/Mb 确认匹配到位：E_s104 = 9,645 vs A_k16 = 9,422（2% 内）。

## 2. 机制（F1_mechanism.tsv，0.5× 中位）

- **每窗观察 landmark 数：sketch 17.3–17.5 vs panel 14.1–21.3，相当**
  （det_frac 0.175 vs 0.185，亦相当）。密度匹配后窗口填充无差异——
  "确定性锚点让低深度窗口保持填充"的说法在本实验体制下**不成立**。
- sketch 在 0.5× 丢 6/17 cell（M1/M10/M13/M18/M23/M25，散布、与生长率
  或覆盖度无关），死因单一：`FracMinHash:no downhill origin: profile
  has no replication gradient`——单 strata 在 0.5× 有 ~1/3 概率拟合层
  建不出梯度，**没有任何融合可救**；panel 的 16 层融合里只要若干层拟合
  成功就能出估计。
- 因此 panel 的低覆盖优势（本体制下）来自**多 strata 融合的冗余**，
  不是 landmark 的确定性。这是对论文 framing 的直接修正，且机制可见
  （每窗 landmark、det_frac、死因三列都支持）。

## 3. 口径警告（必须随结果一起报）

1. **本体制 A 基线 = C8 harness 数字**（A_k16 0.5× r=0.568），不是已提交
   arm A 的 0.913。跨体制比较无效（C8 复核已记录同一现象）；本文件的
   A-vs-E 比较是同体制内部比较，有效。
2. A 臂带 `--windows`（C1 parity）、E 臂不带（armE parity）——各按自家
   流程 parity 跑，不是完全对称 harness。两臂同设 `--windows` 的
   对称版未跑（可作为复核项，预期不改变 §2 机制结论）。
3. sketch 的 r 在幸存 cell 上计算（n=10–11）vs panel 在全 cell 上
   （n=15–16），有幸存者选择效应；但死因（无梯度）表明被丢的是**难**
   cell 而非易 cell，选择效应方向上是保守的。
4. E_s30/s60（密度 3.5×/1.8× 于 panel）为参照臂：0.5× r 0.36（n=5）/
   0.54（n=9），低覆盖下高密度不带来收益（闸门损失扩大为主）；10× 时
   全体打平。密度匹配（而非密度最大化）是正确的问题设定。

## 4. 论文化 consequences（按 HPC_TASKS_SKETCH_MODE.md §5 的预案）

- "deterministic anchors keep windows populated at low depth"（Fig 8
  caption、outline §4、A3 交互项的措辞）**需改写**：匹配密度下窗口填充
  相同；anchors 的剩余优势是湿实验可实现性（route B 不可选）+ 多酶
  融合冗余（低覆盖鲁棒性，机制 = strata 数而非 landmark 性质）。
- 可写论文的新正面结果：密度匹配实验证明 arm E 的低覆盖差距是密度
  伪影；并且定位出 panel 低覆盖鲁棒性的真实机制（融合冗余）。
- F7（混合 landmark）动机增强：两模式失效机制不同（sketch = 拟合层
  无梯度；panel = 单窗 landmark 稀疏），互补有据。
- F2 优先级不变：GC 范围仍是 enzyme panel 作为仪器的边界问题。

## 5. 待办

1. （可选复核）对称 harness：E 臂加 `--windows` 重跑 stats，确认 §2 不变。
2. F2 GC 扫描开工。
3. 若采纳 §4 改写，outline.md Fig 8 caption 与 §4 措辞同步（论文仓侧）。
