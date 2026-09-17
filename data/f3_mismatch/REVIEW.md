# F3 REVIEW — mismatch 容忍下的 landmark 误配风险

日期：2026-09-09。数据：`F3_nearneighbor.tsv`（近邻碰撞普查，108 行 = 18 基因组 ×
2 模式 × tag 长度组）、`F3_misassign.tsv`（实证误配，324 行 = 18 × 2 模式 × 3 错误率
× 3 mm）、`F3_validation.tsv`（Rust↔replica 校验，1782 行）。产物：`fq/`（54 个模拟
FASTQ）、`cells/`（162 个 Rust counts + stats）、`parts/`（逐基因组中间表）、
`scripts/`。索引 db 直接复用 F2 的 72 套 k16 索引，FMH sketch 复用
`bench/F2/sketch_cache/`（k16 匹配档，scale 见 `f2_scales.tsv`）。

作业链（全部 COMPLETED，ExitCode 0）：census 4010749（12s）/ chain 4010750
（17m15s，simreads+counts+replica+aggregate）/ chain 重跑 4010785（9m37s，仅重算
replica+aggregate；原因见口径 6）。调试期 smoke/debug 作业 4010666–4010735 不影响产物。

## 0. 口径（必须随结果一起报）

1. **基因组**：同 F2 的 18 个（GC 25.4–72.0%，`genomes.tsv`）。
2. **普查约定**：对索引 db 中**全部**锚点（含 multi-copy/shared masked 与
   non-chromosomal——它们都可被计数层检索到）枚举同 tag 长度组内的 instance
   对，Hamming 距离 d∈{1,2}（d=0 的相同 tag 对属 multi-copy 家族本身，不是
   mismatch 风险）；跨 tag 长度不对枚举（read tag 长度固定，跨长度不存在混淆通道）。
   酶锚点序列按坐标从参考 FASTA 提取（`contig[position:position+tag_len]` 正链，
   digest.rs 语义），并与 anchors.bin 打包序列逐锚点核对一致（0 不一致）。
   unique = `UNIQUE_IN_GENOME` flag；FMH unique = key 在基因组内仅 1 个 locus
   （armE 口径）。分 uu（双 unique）/ um（恰一 unique）/ mm。
   **普查 ≠ 实证暴露**：普查只覆盖 landmark–landmark 对；FMH 的实证风险主要走
   "rescue 通道"（下详），其 query 不是 landmark，不在普查内。
3. **reads 为自模拟**（`scripts/f3_11_simreads.py`）：F2 模拟器去掉了 PTR 梯度
   （误配不需要），150 bp 单端、0.5×、均匀起点、环化 wrap、50% 反向互补；read 名
   编码真实起点与方向（`@{gid}|{start}|{rc}|{i}`）。错误为 i.i.d. 替换（替换碱基
   均匀 ACGT 抽样，可同碱基），三档：0 / 0.001 / 0.01。**无 indel、无质量偏倚**。
4. **实证误配的实现口径**（任务文档允许的做法，如实报告）：
   - 酶模式：真实 Rust 计数层跑 `--max-mismatch 0/1/2`（162 个 profile，
     `--no-stats`），另用 Python replica 复刻 count.rs 全部计数规则（IUPAC motif
     扫描移植自 enzyme.rs；canonical exact 优先——exact 命中即抑制近邻搜索；
     seed_ranges/best-distance/record-all-best-anchors；keep_multimappers=true），
     利用 read 名中的已知来源坐标把每个被记录的观测分类为 correct /
     wrong_tie（真 locus 与 ≥1 个 Hamming≥1 的邻居同分被记）/ wrong_only
     （邻居严格更近）。**wrong 只计 Hamming≥1 的近邻碰撞**；相同 tag 的
     multi-copy 并列（index 期 masked、EM 层处理）不算 F3 误配。
     replica 与 Rust 在 162 个 cell 上**逐锚点全等**（F3_validation.tsv，
     1782/1782 OK：reads/tag 计数器、mismatch 直方图、multi-locus/enzyme、
     per-anchor counts 全部一致）。
   - FMH 模式：**sketch 计数模式尚不存在于 CLI**（armE 的 count64 只做 exact），
     故为仿真：同一计数规则（exact-canonical 优先、seed 近邻、best-distance）
     作用于 k16 匹配档的 sketch key，k=31 全 read 滑窗（无 motif 门，对应
     Pilea 对每个 read k-mer 计数、per-read 去重的约定；误配率为逐观测口径）。
     即"若 `sk2bgrow index --mode fracminhash` 带 --max-mismatch m 跑，会发生的
     事"。
5. **判定依据**：计数语义读自 `sk2bgrow-core/src/count.rs`（exact-first 抑制、
   best-dist 并列全记、keep_multimappers 默认 true）、`anchor_db.rs`（flags）、
   `tgt.rs`（pack_bases）、`digest.rs`（tag=正链坐标窗口）、`enzyme.rs`（panel）。
   anchors.bin 为 bincode `(Vec<Anchor>, Vec<[u8;12]>)`，26 B/anchor + 12 B tag，
   脚本内解析。
6. **过程事故（不影响产物）**：初版校验脚本的 diff 计数把"表里 0 计数的行"
   误报为 mismatch（`.get(k)` 无默认值的 None≠0）；replica 本身经独立重跑验证
   逐锚全等。修复后 4010785 全量重算，1782/1782 OK。

## 1. 近邻碰撞普查（F3_nearneighbor.tsv）

18 基因组 pooled（853,015 酶锚点 / 855,385 FMH landmark instances，密度匹配 ≤4%）：

| 模式 | d≤1 uu | d≤1 um | d≤1 mm | d≤2 uu | d≤2 um | d≤2 mm | 对/1k landmarks (d≤1) | (d≤2) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 酶板 k16 | 1,328 | 1,934 | 3,224 | 3,191 | 4,438 | 4,972 | 7.60 | 14.77 |
| FMH k16 匹配档 | 13 | 31 | 54 | 31 | 40 | 64 | 0.115 | 0.158 |

- **Syn2b 结论在匹配密度下复现并强化**：FMH 的 landmark–landmark 近邻碰撞比酶板
  低 ~66×（d≤1）/ ~94×（d≤2），d≤2 也几乎为零（0.16 对/1k）。酶板的碰撞集中在
  27 bp 组（9 个酶共享该长度，占锚点 ~57%）。
- **um（unique 贴着 multi-copy 家族）是主要暴露**：占 d≤2 对的 35%（酶）/
  30%（FMH）。per genome 的 d≤2 um 峰值在 GC 33–61%（g09 E. coli 近缘 1,066、
  g05 S. aureus 768、g12 724、g08 452）——与 IS 元件/重复家族丰度一致；两个
  Buchnera（GC 25–26%）几乎为零（基因组仅 0.65 Mb 且几乎无 multi-copy，
  小基因组 confound 与 F2 相同）。
- **GC 非单调性**：d≤2 uu 在两端 GC 反而高（g16/g17 Streptomyces 72%：
  597/360 uu；g14/g15 Mycobacterium 66%：152/194），高 GC 基因组的随机 d=2
  碰撞增多；um 则在中 GC 峰值。两端机制不同。

## 2. 实证误配（F3_misassign.tsv，pooled 18 基因组）

被记录观测中 credited 到 Hamming≥1 错误坐标的比例：

| 错误率 | 酶 mm0 | 酶 mm1 | 酶 mm2 | FMH mm0 | FMH mm1 | FMH mm2 |
|---|---:|---:|---:|---:|---:|---:|
| 0 (e0) | 0 | 0 | 0 | 0 | **9.49e-3** | **2.14e-2** |
| 0.001 | 6.2e-6 | 2.1e-5 | 2.1e-5 | 9.0e-6 | 9.3e-3 | 2.2e-2 |
| 0.01 | 4.1e-5 | 1.0e-4 | 1.4e-4 | 2.2e-5 | 8.2e-3 | 2.0e-2 |

（recorded 总量：酶 268k–328k/格，FMH 272k–350k/格；酶的 wrong 绝对数 worst 为
g09 e01/mm2 的 11/18,099 = 6.1e-4。）

**机制差异（本实验最重要的结构发现）**：

- **酶模式近乎免疫**。count.rs 的 exact-first 规则 + motif 门：无错 read 的 tag
  必有 exact 锚点，近邻搜索被完全抑制；只有"read tag 被错误破坏到无 exact 命中"
  才进入近邻搜索，而 d=1/2 近邻在 9.4k 锚点/Mb 下极稀。0.1% 错误时全 pool
  mm2 仅 7/327,220（2e-5）；1% 错误时 43/316,478（1.4e-4）。
- **FMH 的 mm≥1 风险走 rescue 通道，与测序错误无关**：scale≈100 时 ~99% 的
  read k-mer 本就不在 sketch（hash 阈值过滤），mm≥1 下这些 off-sketch k-mer
  会做近邻搜索，~1%（mm1）/~2%（mm2）撞到一个 sketch key 并被记到**别人的坐标**
  —— e0 时全部 3,274 个 mm1 新增观测（mm1−mm0 的 recorded 增量）无一例外是
  wrong。普查里 FMH landmark–landmandmark 近邻≈0，但 query 空间是 landmark 的
  ~100 倍（全基因组 k-mer），rescue 率 1e-4/观测。这是普查看不到的通道。
- FMH mm0 与酶同样干净（≤2.2e-5，全部来自"错误把 read 变成另一个 sketch key
  的 exact 序列"的罕见事件）。

**敏感性对照（Rust unmatched/motif 窗口，pooled）**：0.1% 错误时 mm0→mm1→mm2
丢失率 2.06%→0.44%→0.42%（mm2 比 mm1 多救回仅 59/328,615）；1% 错误时
19.2%→5.6%→4.4%（mm2 多救 3,693，占 mm0 丢失的 8%）。

## 3. 结论（回答任务的两个决策问题）

1. **`--max-mismatch 2` 是否保持默认？** —— **可以保持；若内存倒逼，降到 1 在
   误配侧是安全的**。误配不是降 mm 的障碍：两种 mm 下误配都 ≤1.4e-4（1% 错误），
   0.1% 错误时 mm1 与 mm2 误配相同（7/327k）。mm2 的真实代价只在高错误率 reads
   上买回 ~8% 的可恢复观测；在 0.1% 错误下 mm2 相对 mm1 几乎无收益。结合 A1 的
   索引内存（mm0→mm1 实测 174→356 B/anchor，seed slot 随 mm+1 增长），**若
   C6 的 GTDB 预算（752 GB @16 酶）是硬约束，默认降 mm=1 误配侧无代价、灵敏度
   损失可忽略（0.4% 窗口丢失）；reads 质量差（如考古/ FFPE）时再切回 mm=2**。
   这是工程-预算决策，F3 只负责把误配风险从决策里摘掉。
2. **sketch 模式能否安全用更高 tolerance？** —— **不能；sketch 模式应锁
   mm=0**。mm1 就有 ~1%、mm2 ~2% 的观测被记到错误坐标（酶模式的 ~100–1000×），
   且与测序错误率无关（rescue 通道）。若未来给 `--mode fracminhash` 暴露
   `--max-mismatch`，应默认 0 并警告；或者只在 query 同时满足"参考该位有
   landmark"时才允许近邻救援（当前计数语义天然满足前者于酶模式，FMH 无此门）。
   换言之：Syn2b"FMH 近邻碰撞 0.00%"的说法在 landmark–landmark 意义上成立
   （本普查确认），但**不能**推出"sketch 模式可以用高 tolerance"——恰恰相反。

## 4. 未决与警告

- FMH 实证为规则仿真（sketch 计数模式不存在于 CLI）；酶模式为真实 Rust 计数
  + 全等校验的 replica 归因。两者口径已尽量对齐（同计数规则、同 best-dist
  语义），但 FMH 侧没有独立实现可校验。
- 错误模型仅替换；indel/质量偏倚未测。1% 档是压力测试，不代表常规数据。
- 0.5× 深度、单端 150 bp；深度不改变误配机制（每观测概率），未扫深度。
- 普查含 non-chromosomal 锚点（可检索即暴露）；若只关心染色体建模，um/uu
  数会略降，不改结论。
- Buchnera 的低碰撞部分是小基因组/低重复 confound，不应读成"低 GC 安全"。
- F3 不测灵敏度收益本身（那是 A1/C6 的预算问题）；§2 的 unmatched 表只是
  决策参照。
