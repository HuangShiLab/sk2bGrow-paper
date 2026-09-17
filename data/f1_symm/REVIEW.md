# f1_symm REVIEW — F1 对称 harness 重跑（--windows flag 判定）

日期：2026-09-12。数据：`F1_symm_results.tsv`（2 flag 配置 × 3 arms × 5
深度）、`F1_symm_all_cells.tsv`（逐 cell）、`out/<flag_config>/<arm>_<medium>_<cov>x/`。
作业：sketch250 = 4025242（amd，~20 min）/ stats = 4025243（amd，~13 min，
--dependency=afterok:4025242）。两者均 0 FAIL 0 MISSING 完成。

## 0. 设计（与原始 F1 的差异点）

- reads/counts 完全复用原始 F1 同一批：`$F1/cells/k16/`（85 cells）、
  `$F1/sketch/s104/`（85 cells）；E_s250 的 counts 为本目录新建
  （`sketch/s250/`，armE_counts.py，cache 自建于本目录，未动 F1 文件）。
- db 复用 `$BASE/bench/C1/db`（A_k16 与两个 E 臂同原始 f1_13）。
- 两种 flag 配置：**with_windows**（A 传真实 `windows.tsv`，E 传自身
  counts.tsv 作占位 —— 见 §1）与 **without_windows**（均不传）。
- arms：A_k16、E_s104（与 panel 密度匹配，9645 vs 9422 landmarks/Mb）、
  E_s250（新，3998 landmarks/Mb，~0.42× panel 密度）。

## 1. 关键技术发现：`--windows` 在当前代码里是 no-op

`sk2bgrow.cli profile` 解析 `--windows`（cli.py:156）但 `_profile` 从不读
`args.windows`；真正起作用的开关是 `--use-rust-windows`（cli.py:47,162，
经 `use_precomputed_windows` 传入 ztp.window_rates）。Rust 侧
`run_python_stats`（profile.rs:368）同样只传 `--windows`、不传
`--use-rust-windows`。原始 F1 的"A 带 --windows、E 不带"在 effect 上是
**对称的**，不是 0.568 vs 0.913 差异的来源。

实证：本重跑 255 对 with/without 输出目录的 `output.tsv` 逐字节
`cmp`，**255 identical / 0 differ**；逐 cell windows.rates.tsv 同
（数字来自同一 count 表、同一 db、同一 fitting 路径）。

## 2. 三个判据的回答

### (a) A_k16@0.5× 是否复现 0.568？—— 是，精确复现

| flag_config | A_k16 @0.5× r（n） |
|---|---|
| with_windows | **0.5681（16）** |
| without_windows | **0.5681（16）** |
| 原始 F1_results.tsv | 0.5681（16） |

逐深度全部与原始 F1 一致（1× 0.9060、2× 0.9349、5× 0.9797、10×
0.9854；E_s104 亦完全一致）。→ **0.568 vs C8 的 0.913 是 harness 固有
差异，不是 --windows flag artifact**。真实来源：C8 头条数字是 GLM 拟合
方法（C8 CLOSEOUT：0.5× glm r=0.911 vs v_shape 0.568），F1 用默认
v_shape。计时/比较口径提醒不变：跨 harness 数字不可比。

### (b) 0.5× 存活细胞上 A-vs-E 排序是否翻转？—— 不翻转

| flag_config | 臂 | r（n=存活） | 同 cell 集合 A_k16 r |
|---|---|---|---|
| with / without | A_k16 | 0.5681（16，全集） | — |
| with / without | E_s104 | **0.8183（10）** | 0.6905（同 10 cell） |
| with / without | E_s250 | 0.5527（13） | 0.4828（同 13 cell） |

- E_s104 死 cell 集合 {M1,M10,M13,M18,M23,M25} 与原始 F1 完全一致；
  死因与原始一致（拟合层建不出梯度，见 stats.log）。
- 密度匹配的 E_s104 在幸存 cell 上仍显著高于 panel——**即使把 A 限制
  到同一批 10 个幸存 cell（A 也只有 0.6905），E_s104（0.8183）仍领先**。
  对称 flags 不改变排序。
- 新臂 E_s250（低密度，3998 landmarks/Mb）：全集 r=0.5527（13）≈ 打平
  略低于 A_k16 的 0.5681（16）；同 cell 集合上 E_s250 0.5527 > A 0.4828。
  即"sketch 在存活细胞上反超 panel"的现象在密度匹配（s104）时成立，
  在 0.42× 密度（s250）时基本消失——与 F1"密度匹配是关键设定"的结论
  同向。

### (c) ≥1× "两 landmark 源不可区分"是否不变？—— 不变

| 深度 | A_k16 r（n） | E_s104 r（n） | E_s250 r（n） | 判定（两配置数字相同） |
|---|---|---|---|---|
| 1× | 0.9060（16） | 0.9089（16） | 0.8391（16） | A≈E_s104 打平 |
| 2× | 0.9349（16） | 0.9662（16） | 0.9478（16） | 打平（Δ≤0.03） |
| 5× | 0.9797（16） | 0.9863（16） | 0.9819（16） | 打平 |
| 10× | 0.9854（16） | 0.9812（16） | 0.9799（16） | 打平 |

与原始 F1 的 ≥1× 结论完全一致；E_s250 在 1× 略低（0.8391）但同量级，
不推翻"不可区分"。

## 3. 对 F1 REVIEW 结论的影响

- F1 REVIEW §3.2 的口径警告（"A 带 --windows、E 不带，非完全对称"）**降
  级为 cosmetic**：当前代码 --windows 无 effect（§1），原始 F1 的
  A-vs-E 比较本身就是有效的对称比较。
- F1 §2 机制结论（窗口填充相当、sketch 死因=拟合层无梯度、panel 冗余
  来自多 strata 融合）不受影响：死 cell 集合、r 值逐位复现。
- 新增可写点：密度降低 2.4×（s104→s250）使 0.5× 幸存 cell 的 sketch
  优势从 +0.13 r（同 cell）缩到 +0.07 r 且全集打平——"sketch 低覆盖
  优势"对密度敏感，密度匹配仍是正确设定。

## 4. 产物与复现

- 脚本：`f1s_01_sketch250.sh` / `f1s_02_stats.sh` / `f1s_03_aggregate.py`
  （本目录）；日志 `$BASE/logs/f1_symm/`。
- 复用未重建：`$F1/cells/k16`、`$F1/sketch/s104`、`$BASE/bench/C1/db`。
- 聚合：`sbatch 4025242 → 4025243` 后跑 f1s_03_aggregate.py；任何缺失
  cell 在结果表中以 n 减少体现，可补跑 stats（skip-if-exists 幂等）后重
  聚合。
