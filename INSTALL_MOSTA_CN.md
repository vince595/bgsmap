# bgsmap v1.0.0 在 MOSTA E16.5 数据上的安装与使用

适用于已有 gsMap 1.73.7 和 MOSTA GSS 文件的环境。可复用已生成的 GSS。

## 1. 在服务器上安装新 WHL

在保存 `bgsmap-1.0.0-py3-none-any.whl` 的目录中执行以下命令。

```bash
conda activate gsmap
python -m pip install './bgsmap-1.0.0-py3-none-any.whl[all]'

python -c 'import bgsmap; print(bgsmap.__version__, bgsmap.__file__)'
bgsmap --version
```

预期版本为 `1.0.0`。旧项目的分发名为 `gsmap-burden`，新项目的分发名、导入名和 CLI 均为 `bgsmap`。

也可从源码 ZIP 安装。

```bash
unzip bgsmap-v1.0.0-source.zip
cd bgsmap-1.0.0
python -m pip install '.[all]'
```

两种方式**二选一**；推荐 WHL。

## 2. 先执行无大矩阵加载的 spot 匹配检查

```bash
BGS=/path/to/project
ROOT=/path/to/dataset
WORK="$BGS/Height_MOSTA_E16/GWAS/gsmap_work"
SAMPLE=E16.5_E1S1.MOSTA
ST="$ROOT/data/gsmap/mouse-e16-mosta/ST/E16.5_E1S1.MOSTA.h5ad"
GSS="$WORK/$SAMPLE/latent_to_gene/${SAMPLE}_gene_marker_score.feather"
BURDEN="$BGS/Height_MOSTA_E16/Burden/input/Height_pLoF_QC.tsv"

bgsmap inspect \
    --gss "$GSS" \
    --h5ad "$ST" \
    --spatial-key spatial \
    --annotation annotation \
    --check-spots
```

预期 `"matched": 121767`。新增的 `inspect` 对 Feather 使用 `pyarrow` 的 schema-only 读取，不会将 121,767 个 spot 列全部转换成 pandas 矩阵。

## 3. 测试运行（解析 P，0 次置换）

```bash
OUT="$BGS/Height_MOSTA_E16/Burden/Height_pLoF_smoke"

bgsmap run \
    --gss "$GSS" \
    --gss-gene-col HUMAN_GENE_SYM \
    --h5ad "$ST" \
    --spatial-key spatial \
    --annotation annotation \
    --burden "$BURDEN" \
    --gene-col gene_symbol \
    --p-col Pvalue_Burden \
    --trait-name Height_pLoF \
    --score-mode logp \
    --p-cap 30 \
    --permutations 0 \
    --chunk-size 32 \
    --seed 42 \
    --out "$OUT"
```

如果成功，检查 `run_metadata.json` 中的 `n_spots` 和 `n_genes`。正式主 P 值不依赖 permutation 次数；只有需要额外校准时才使用 `--permutations 999` 等参数。121,767 个 spot 与大量基因的密集矩阵计算仍可能占用较高内存。

## 4. 正式运行（解析 P，0 次置换；按需增加校准置换）

```bash
OUT="$BGS/Height_MOSTA_E16/Burden/Height_pLoF"

bgsmap run \
    --gss "$GSS" \
    --gss-gene-col HUMAN_GENE_SYM \
    --h5ad "$ST" \
    --spatial-key spatial \
    --annotation annotation \
    --burden "$BURDEN" \
    --gene-col gene_symbol \
    --p-col Pvalue_Burden \
    --trait-name Height_pLoF \
    --score-mode logp \
    --p-cap 30 \
    --permutations 0 \
    --chunk-size 32 \
    --seed 42 \
    --out "$OUT"
```

**注意：** v0.2.0 已重新设计 P 值与全部主图：解析 P 为主，经验 P 仅保存于可选诊断文件（见 `docs/PVALUES_v0.2_CN.md`）。`--match-cols` 可用于条件基因可交换性建模。若使用 `--coords`，仍不自动启用 AnnData 回退。

## 5. 如果仍报 `Only 0 GSS spots...`

请提供：

```bash
python -c 'import bgsmap; print(bgsmap.__version__, bgsmap.__file__)'
bgsmap inspect --gss "$GSS" --h5ad "$ST" --annotation annotation --check-spots
```

新异常会包括 GSS/H5AD 读取后的前五个 ID，可精确定位是哪个适配器产生了差异。不要在没有逐项确认 spot ID 对应关系时按行顺序强行拼接。
