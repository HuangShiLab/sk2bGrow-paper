# B/C 两臂提取规则保真核查（2026-09-12）

## 1. arm C 提取规则 —— 无系统偏差

- `extract_armC.py` 取每个 pilea output.tsv 的**首行**（`next(csv.DictReader)`），字段映射 `log2(PTR)`→log2ptr、`coverage`→est_cov；代码不做基因组名过滤，但 pileadb 是单基因组库（`genomes.tab` 仅 `Escherichia_coli_K12`，k=31/s=250/w=25000），且已核验 **全部 630 个 output.tsv（45 runs×7 深度×2 模式）零多行文件**——每个文件要么仅表头（gate 拒绝），要么表头+唯一一行 → 首行提取在结构上无歧义，不存在"取错质粒/重叠群/他基因组行"的可能。
- 抽查：M6@2× relaxed 原始行 = E. coli K12、log2(PTR)=1.4398、coverage=1.9033；M10@5× relaxed = 1.2734/3.6328，与 `regen_armC_cells.tsv` 逐位一致；default 模式两格均为表头（gate 拒绝）→ 估计留空，与 Mac 备份表示法一致。
- 版本：Pilea **1.3.8**（`pilea --version` 与 pip show 一致，同论文声明）；参考序列 `refs/Escherichia_coli_K12.fna` 头行 = NC_000913.3（E. coli K-12 substr. MG1655，即 GCF_000005845.2）。

## 2. arm B（sorted）保真 —— 确认为 RANSAC sorted-rank 回归

- M10@2× 完整输出行：method=**sorted_ransac**、n_windows=539、log2(PTR)=1.016815635，与 `regen_armB_cells.tsv` 逐位一致；**85/85 输出 method 列全部为 sorted_ransac**，无静默回退。
- 实现位置：`src/python/sk2bgrow/fit.py:115` `fit_sorted_ransac`（文件头注释 5–6 行、docstring 122–124 行自述 "iRep/Pilea-style sorted regression with RANSAC"：窗口 log2 rate 排序、对 rank 回归、斜率×全距）。
- 与 arm A 同窗口表：同 counts、同 `--windows`（C1 windows.tsv），M10@2× 的 A/B n_windows 均为 539，仅 fit 方法不同（v_shape vs sorted_ransac）。
- 交叉验证：M10@2× B=1.0168 vs A=0.9960（真值 1.15），差 0.021，量级合理。

## 3. 方向性检验 —— B 为全局负向平移，C 为温和散点

- **arm B（regen−Mac，n=85）**：分深度均值全负——0.5× −0.258、1× −0.666、2× −0.549、5× −0.278、10× −0.151；mean|d|=0.387；|d|>0.5 共 28 格，>1.0 仅 1 格（M17@1× −1.008，刚过线）。形态是"全网格一致的下移"，与提取 bug 应有的"稀疏剧烈离群"不符 → 判定为**估计器代码漂移**（Mac 8/25 实例与现行工作树之间 fit 层变动），非提取错误。
- **C_relaxed**：mean|d|=0.096，**无 |d|>0.5 格**，0.5× 两侧逐位相同（全 0）；2× r 下滑（Mac 0.947→regen 0.848）源自散点漂移（最大 M13@2× −0.476、M3@2× −0.351），无单格异常。

## 4. 对交互变化归因的影响判定

B、C 两臂提取规则均排除系统性偏差：C 是逐字节忠实的搬运；B 是真实的 sorted 估计且与 A 同窗口。2×2 交互检验的变化（旧 1×/2×/10× 正显著 → 新全不显著或 2× 负）应维持归因于**实例差异**——具体为统计层代码漂移（B 臂系统性下移 + C 温和散点 + 此前已证的 A 臂逐 cell 漂移），其中 B 的下移恰好同时削弱 estimator 主效应与交互项，属于漂移的真实组成部分而非 artifacts。建议：Table 2 报告新实例数字，并在方法部分记录估计器代码版本（git describe）与运行日期；若需追查到具体 commit，需找回 Mac 8/25 的工作树快照（本机已无）。
