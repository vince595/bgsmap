# bgsmap 的连续 P 值、结果保存与绘图

## 设计原则

- **主分析不再使用 `p_empirical`**：所有绘图、排序和 `region_results.tsv` / `spot_results.tsv.gz` 的主 P 值使用 `p_analytic` / `p_competitive_analytic`。
- **不受 1/(permutations+1) 的蒙特卡洛分辨率限制**：采用分层基因标签置换零假设下精确的二阶矩，并以正态近似计算单侧尾概率。
- **极小 P 不受普通浮点存储下限限制**：保存 `logp_analytic` / `logp_competitive_analytic`，它们是 `-log10(P)` 的数值，而不是对存为 0 的 `p` 再取 log。
- **原有单表型内 BH 校正不扩大范围**：输出 `q_bh_analytic` / `q_bh_competitive_analytic`。额外给出单表型内 Holm FWER 作为参考；不实施跨 108 表型联合校正。
- **permutation 为可选校准诊断**：`--permutations 0` 为默认。使用 999/9999 次时，旧经验 P、maxT-FWER 等仅保存到 `*_permutation_diagnostics.tsv*`，不会影响图或主结果文件。

## 统计原理（重要）

对于基因 `g`、spot `s`，将该基因 across spots 的 GSS 标准化为 `X_gs`；从 burden P 计算强度 `w_g`，并按指定匹配分层 `h` 进行去均值和全局缩放。spot 统计量：

`T_s = (1/G) sum_g X_gs w_g`。

对于每个 stratum 中固定的 `X` 和固定的 `w` 值，将 `w` 随机分配给该层的基因，其条件零分布具有精确方差：

`Var(T_s) = G^{-2} sum_h [ sum_{g in h}(w_g - mean_h w)^2 / (n_h - 1) ] * sum_{g in h}(X_gs - mean_h X_s)^2`。

于是令 `Z_s = T_s / sqrt(Var(T_s))`，用 `P_analytic = Pr(N(0,1) >= Z_s)` 作为**连续的单侧正态近似 P 值**。

同样分别用域平均 GSS 和域内减域外平均 GSS 计算域 `T_d` 和竞争性域 `T_comp,d`，按同一解析方差公式求 `Z` 和 `P`。

**局限**：方差公式在分层可交换假设下是精确的，但**正态尾部近似不是精确 P 值**，并不保证非常小的 P 值一定校准；特别是 burden 权重由少数基因主导、基因间相关性强、功效差异未匹配、或 GSS 具有极端值时。`run_metadata.json` 的 `effective_gene_weight_count` 和 `normal_approximation_warning` 提供基本预警；需要在实际论文的 global null / competitive null 模拟中检验 P 值校准。不能将连续 P 值等同于已被验证的极端尾部概率。

## 主要输出

| 结果文件 | 字段 | 用途 |
|---|---|---|
| `spot_results.tsv.gz` | `score`, `z_analytic`, `p_analytic`, `logp_analytic`, `q_bh_analytic`, `logq_bh_analytic`, `p_holm_analytic` | spot 统计与空间图着色 |
| `region_results.tsv` | `mean_score`, `p_analytic`, `logp_analytic`, `q_bh_analytic` | domain mean 富集 |
| `region_results.tsv` | `competitive_score`, `z_competitive_analytic`, `p_competitive_analytic`, `logp_competitive_analytic`, `q_bh_competitive_analytic`, `p_holm_competitive_analytic` | **竞争性域富集的主结果和域图** |
| `spot_permutation_diagnostics.tsv.gz` | `p_empirical`, `q_bh_empirical`, `p_fwer_permutation` | 仅在 `--permutations > 0` 时存在，校准和敏感性检查 |
| `region_permutation_diagnostics.tsv` | `p_empirical`, `p_competitive_empirical`, `p_fwer_permutation`, `p_competitive_fwer_permutation` 等 | 仅在 `--permutations > 0` 时存在，校准和敏感性检查 |
| `run_metadata.json` | `primary_p_column`, `primary_competitive_p_column`, `effective_gene_weight_count`, `primary_method` | 方法及近似可靠性声明 |

注意：`p_analytic` 可能因 IEEE 754 数值下溢存为 `0`，并不意味着精确 P=0。图形、极强信号排名和可视化必须使用 `logp_analytic`，不要对存储的 `p_analytic` 直接求 `-log10`。

## 绘图变更

- `plots/spatial_enrichment.html` / `.pdf` 左图：`logp_analytic`，标题 **Analytic enrichment (-log10 p)**。右图：`score`。
- `plots/region_enrichment.html` / `.pdf`：**`logp_competitive_analytic`**，不再计算 `-log10(p_empirical)` / `-log10(p_competitive)`。
- `plots/report.html`：按 `p_competitive_analytic` 排序，表中展示解析 P、Z、BH 和 Holm；不把 permutation 经验值混进主表。
- `plots/gene_contributions.pdf`：仍为描述性贡献图，主域由 `p_competitive_analytic` 决定。

## 服务器安装（清华镜像）

```bash
# PyPI 完成发布且镜像同步后
python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple 'bgsmap[all]'
python -m bgsmap --version

# 纯解析 P，无需 permutation
# 将 GSS、BURDEN、H5AD、ANNOTATION、OUT 替换为实际输入
# bgsmap run --gss "$GSS" --h5ad "$H5AD" --annotation "$ANNOTATION" \
#   --burden "$BURDEN" --trait-name Height --permutations 0 --out "$OUT"

# 如需独立诊断，另加 --permutations 999（不会限制主要输出 P）
# 已有解析 P 值结果时，可直接重新生成图
# bgsmap plot --out "$OUT" --trait-name Height
```

**旧结果迁移**。v0.1.x 生成的 `spot_results.tsv.gz` 和 `region_results.tsv` 只有经验 P。请使用同样的输入、基因映射和参数重新运行 `bgsmap run` 以获得解析结果。`bgsmap plot` 会在缺少解析 P 值时提示重新运行。
