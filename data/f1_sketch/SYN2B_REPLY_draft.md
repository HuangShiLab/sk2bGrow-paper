# F6 — Read-level FracMinHash sketch: answer to Syn2b §4 (draft)

Syn2b（LANDMARK_COMPARISON）§4 列的未决问题："Behaviour on reads rather than
assemblies. FracMinHash needs an assembly. Whether a read-level sketch
reproduces these results is untested."

sk2bGrow 的 F1/arm E 就是这个实验，且比"reproduce"更强：不是把 sketch
用在组装级丰度上，而是

1. **landmark 定义在组装参考上，计数来自 150bp reads**（0.5–10× 深度梯度）；
2. 经过一个**坐标依赖**的估计器（log2 窗口速率对离 ori 的环形距离回归 +
   原点搜索 + 融合），比纯丰度度量对 landmark 质量更敏感；
3. 在真实数据上有独立真值（Zheng et al. 16 培养条件实测生长率）。

结果（密度匹配后，E. coli K12，landmarks/Mb 9,645 sketch vs 9,422 酶板，
差 2% 以内）：

| 深度 | sketch r（n） | 酶板 r（n） |
|---|---|---|
| 0.5× | 0.82–0.84（10–11） | 0.50–0.57（15–16） |
| 1× | 0.85–0.92（16–17） | 0.89–0.91（15–16） |
| 2× | 0.92–0.97 | 0.93–0.98 |
| 5×/10× | ≈0.98/≈0.98 | ≈0.98/≈0.98 |

GC 泛化（F2，18 基因组 spanning 25–72% GC，planted-gradient 模拟）：
匹配密度下 sketch 与酶板在全 GC 范围打平（25–35% 低 GC 带 0.5× sketch
甚至略好）；仪器边界是 GC≲30% + 0.5×，与 landmark 类型无关。

对 Syn2b 问题的直接回答：**read-level sketch 成立**。剩余的两点限定
（诚实声明）：(a) 0.5× 时单 strata sketch 有 ~1/3 的 cell 在拟合层建不出
复制梯度（no downhill origin）被闸门拒绝，酶板靠多 strata 融合冗余救回
——read-level sketch 的精度要乘以存活率来看；(b) 我们的 read 模拟是
无测序错误的，Syn2b 的 adjacency 问题（read-level sketch 是否保住
landmark 邻接信息）我们未测。

证据出处：bench/F1/REVIEW.md、bench/F2/REVIEW.md（sk2bGrow 仓库）。
