# bgsmap 中文使用手册

bgsmap 联合原版 gsMap 产生的基因特异性评分（GSS）和基因级稀有变异 burden 汇总统计，计算 spot 和空间域富集，输出 TSV、交互式 HTML 和矢量 PDF。

主结果和绘图使用解析 P 值。默认 `--permutations 0`，也可按需开启独立的经验 P 值诊断。解析 P 值基于正态尾部近似，研究中使用极小 P 值前应评估校准情况。[P 值与绘图说明](docs/PVALUES_v0.2_CN.md)

1.3.0 将通用分块引擎接入 `bgsmap run`。解析模式自动启用单表型内部 spot 并行，
支持 Feather V2、Parquet、TSV、CSV，按内存预算调整块大小，并记录各阶段耗时。
可复用同一切片的磁盘缓存。详见 [大数据性能与内存指南](docs/PERFORMANCE_CN.md)。

## 1. 安装

要求 Python 3.10 或更高版本。使用以下任一命令安装或升级，清华源需等待镜像同步。

```bash
python -m pip install --upgrade 'bgsmap[all]'
python -m pip install --upgrade -i https://pypi.tuna.tsinghua.edu.cn/simple 'bgsmap[all]'
```

也可安装发布 wheel，或进入解压后的源码目录安装。

```bash
python -m pip install ./bgsmap-1.3.0-py3-none-any.whl
# 或在解压后的源码根目录执行
python -m pip install .
```

读取 Feather 和 H5AD 时，请从 PyPI 安装 `'bgsmap[all]'`，或从源码安装 `'.[all]'`。`feather` 扩展也支持读取 Parquet。AnnData spot ID 回退读取需要额外安装 `anndata`。从 H5AD 训练 GNN 和计算 GSS 的 `quick` 命令需要另装 gsMap。

```bash
bgsmap --version
python -m bgsmap --help
```

## 2. 直接用你的原有 GSS 和 FinnGen burden（建议先跑这个）

gsMap 官方旧版的 GSS 一般位于：

```text
WORKDIR/SAMPLE/latent_to_gene/SAMPLE_gene_marker_score.feather
```

`Feather` 第一列应为 `HUMAN_GENE_SYM`，其他列为 spot barcode；原始 .h5ad 中 `obs_names` 须一致。

```bash
GSS='/真实路径/VisiumCoronal_gene_marker_score.feather'
ST='/真实路径/V1_Adult_Mouse_Brain_Coronal_Section.h5ad'
BURDEN='/真实路径/FinnGen_R12_G6_MIGRAINE_Mask1.0.01_ADD.tsv.gz'

# 先看数据有哪些列，避免注释、P列写错
bgsmap inspect --h5ad "$ST" --burden "$BURDEN"

bgsmap run \
  --gss "$GSS" \
  --h5ad "$ST" \
  --annotation YOUR_DOMAIN_COLUMN \
  --burden "$BURDEN" \
  --gene-col YOUR_GENE_COLUMN --p-col YOUR_P_COLUMN \
  --trait-name FinnGen_G6_MIGRAINE \
  --permutations 0 --seed 42 \
  --out results/FinnGen_G6_MIGRAINE
```

把 `YOUR_DOMAIN_COLUMN`、`YOUR_GENE_COLUMN`、`YOUR_P_COLUMN` 换成检查到的实际列名。若原来的 H5AD 还没有 domain 注释，请自行导出 4 列 `spot,x,y,domain` 的 TSV 文件，以 `--coords` 代替 `--h5ad`。建议提供空间域标签以解释 competitive null。

**关键点：鼠脑 GSS 与人类 FinnGen burden 必须经过严格的基因匹配。** 原版 gsMap 若已经根据 homolog 转换到人类符号，直接匹配即可；否则 `--gene-map` 提供一对一鼠–人基因映射文件，含 `st_gene` 和 `burden_gene` 两列。

如果 FinnGen 表格有多个 mask / trait，指定 `--mask-col` + `--mask`，或 `--trait-col` + `--trait`；切勿简单地按基因直接取最小 P 却不调整多重检验。

## 3. 显式匹配每个基因的 burden 检验能力（推荐用于论文）

burden P 的分布可能和基因长度、可检验变异个数、累计 MAC、测序质量有关。**默认不匹配时只能视为探索性分析。** 可提供基因层面检验机会/功效相关协变量，执行分层置换：

```bash
bgsmap run \
  --gss "$GSS" --h5ad "$ST" --annotation YOUR_DOMAIN_COLUMN \
  --burden "$BURDEN" --gene-col YOUR_GENE_COLUMN --p-col YOUR_P_COLUMN \
  --covariates /真实路径/gene_power_covariates.tsv \
  --match-cols gene_length n_tested_variants \
  --match-bins 4 --permutations 0 \
  --out results/FinnGen_matched
```

协变量文件 `gene,gene_length,n_tested_variants`。若这些列已在 burden 表里，省略 `--covariates`。需要自己确认这些列的定义与 mask 匹配，且每个组合分层至少 12 个基因。gene–gene 相关性、样本量和检测能力异质性仍可能导致校准问题，需独立模拟评价。

## 4. 从原始 H5AD 开始跑（使用原版 gsMap）

```bash
pip install gsMap

bgsmap quick \
  --h5ad "$ST" \
  --workdir work/gsmap_gss --sample-name VisiumCoronal \
  --data-layer counts --annotation YOUR_DOMAIN_COLUMN \
  --homolog-file /真实路径/mouse_human_homolog.tsv \
  --burden "$BURDEN" \
  --gene-col YOUR_GENE_COLUMN --p-col YOUR_P_COLUMN \
  --out results/VisiumCoronal --permutations 0
```

前两步调用原版 `gsMap.find_latent_representation` 和 `gsMap.latent_to_gene`；之后执行本包的 burden 富集与可视化。若 H5AD 中已有原版 GNN embedding，可尝试 `--latent-key latent_GVAE` 跳过 GNN 训练。`--data-layer` 必须指向实际 raw count 层或 `X`。**由于无法访问你的服务器环境，此桥接尚未在你的 gsMap 安装版本和真实 H5AD 上执行验收。**

## 5. 输出说明

```text
results/VisiumCoronal/
├── aligned_burden.tsv.gz             # 匹配后的 burden P，p_original 保留输入数值
├── spot_results.tsv.gz               # 逐 spot 解析 P、logP、BH q 与 Holm
├── region_results.tsv               # 域平均及 competitive 解析 P、Z、logP、BH / Holm
├── gene_contributions.tsv.gz         # 贡献基因，仅为描述性结果
├── run_metadata.json                 # 输入、匹配数、随机种子、统计假设
├── spot_permutation_diagnostics.tsv.gz  # 仅设置 permutations>0 时生成
├── region_permutation_diagnostics.tsv   # 仅设置 permutations>0 时生成
└── plots/
    ├── report.html                   # HTML 总报告
    ├── spatial_enrichment.html       # 可缩放、悬停的空转富集地图
    ├── spatial_enrichment.pdf        # 可编辑矢量点位 PDF
    ├── region_enrichment.html
    ├── region_enrichment.pdf
    └── gene_contributions.pdf
```

`spatial_enrichment.html` 为离线自包含 HTML；`report.html` 使用相邻 HTML 作为 iframe，复制报告时请保留整个 `plots/` 目录。PDF 使用矢量点和可编辑字体；在 Mac Adobe Illustrator 中可打开，字体取决于本机安装情况。

三类统计量：**spot 富集**指此 spot 的 GSS 和真实 burden 强度是否高于基因标签置换；**domain mean 富集**指整个域的平均评分是否异常偏高；**competitive domain 富集**指该域相对于其他 spots 是否表现出更高的 burden 信号。最后一个更适合检验你之前 C（competitive null）的目的。主结果使用解析 P 值 `p_analytic`/`p_competitive_analytic`，`q_bh_analytic` 为单表型内 BH，`p_holm_analytic` 为 Holm 校正；只有可选诊断文件保留原 maxT 与经验 P。

`logp_analytic` 不受 permutation 次数的分辨率限制。`--permutations` 用于额外校准诊断。

## 6. 用模拟数据自测

```bash
# 生成纯模拟数据、跑 HTML/PDF，检查能否在你服务器上运行
bgsmap demo --out demo_data
bgsmap run --gss demo_data/GSS.tsv.gz \
  --burden demo_data/burden.tsv.gz --coords demo_data/spatial_coords.tsv \
  --permutations 0 --out demo_results

python -m pip install -e '.[all,dev]'
python -m pytest -q
python -m build
```

命令默认将带时间戳的运行日志写到 stderr，`run` 和 `inspect` 的 JSON 输出仍走 stdout。需要分别保存时，在 `run` 命令后加 `> summary.json 2> run.log`。

论文中使用 gsMap 计算的 GSS 时，请引用原版论文。[Song et al., *Nature* 641, 932–941 (2025).](https://www.nature.com/articles/s41586-025-08757-x)

### 发布新版本

首次发布时，在 PyPI 账户的 Publishing 页面添加 GitHub Trusted Publisher。项目名填 `bgsmap`，Owner 填 `vince595`，Repository 填 `bgsmap`，Workflow name 填 `publish.yml`，Environment name 填 `pypi`。[PyPI 首次发布说明](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/)

将版本号更新合并到 `main` 后会自动运行 **Publish to PyPI**，也可在 GitHub Actions 手动运行。工作流先执行测试和包检查，再发布到 PyPI，并创建同版本 GitHub Release，附带 wheel 和源码包。

## 7. 服务器兼容性提示

- 建议设置 `export OPENBLAS_NUM_THREADS=4`（按节点实际 CPU 调整），避免矩阵置换时多线程过度占用。`--chunk-size 64` 决定每次处理的置换批量，内存紧张可调为 16。
- 如果你的 `pyarrow` 报 `GLIBCXX_*` 动态库错误，可先在干净的 conda 环境中安装兼容的 `pyarrow` 和 `libstdcxx-ng`。软件会尝试 `pyarrow.ipc.open_file` 作为 Feather V2 读取后备，但**不能保证解决底层 C++ 库缺失**；也可在能读取 Feather 的环境先转换为 `.tsv.gz`。
- 不要把模拟过程人为注入 Z-shift 后重新计算的 P 值描述为“原始 P 值未变化”。原始 P 仅用于 baseline；最终被注入的基因 P 已经变化。本包在真实分析中不会改写输入 P，只在内部构造基因标签置换的 null。


## v0.1.1 历史修复：MOSTA E16.5（121,767 bin）H5AD/GSS spot 匹配

旧版 `gsmap-burden v0.1.0` 使用低层 `h5py` 读取 H5AD spot ID，
并在匹配为 0 时直接失败。新版本在该读取器匹配不到 spot 时，会
尝试通过 `anndata.read_h5ad(..., backed="r")` 恢复标准 `obs_names`，
不会修改原始文件、GSS 或坐标；未恢复时显示双方前五个 ID 供排查。

安装（保持当前 gsMap 环境）：

```bash
conda activate gsmap
python -m pip install --no-deps --force-reinstall bgsmap-1.3.0-py3-none-any.whl
python -c 'import bgsmap; print(bgsmap.__version__)'
```

诊断（只读取 Feather 文件的 Arrow schema，不读取大矩阵）：

```bash
bgsmap inspect \
  --gss /path/to/E16.5_E1S1.MOSTA_gene_marker_score.feather \
  --h5ad /path/to/E16.5_E1S1.MOSTA.h5ad \
  --annotation annotation \
  --check-spots
```

预期 `spot_alignment.matched` 为 `121767`。如果与用户用 AnnData 检查
得到的数量不一致，请将 inspect 的完整输出用于进一步诊断；**不要将
不同文件的 spot 强制按行顺序匹配**。

这次是 IO/可诊断性补丁，**没有修改**基因置换检验公式或空间统计模型。
对于大规模 MOSTA 数据，先用 `--permutations 0` 进行连通性
测试；是否需要额外做 `999` 次 permutation 校准，可再单独决定。
