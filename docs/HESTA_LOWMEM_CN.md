# bgsmap 1.2.0 · HESTA 大切片运行记录

当前通用版本为 1.3.0。安装、自动分块、内存预算及跨数据集用法见
[大数据性能与内存指南](PERFORMANCE_CN.md)。本文保留 1.2.0 的测量记录和集群示例。

1.2.0 将低内存分析集成到 `bgsmap`。使用 `bgsmap-lowmem run` 或
`bgsmap run --engine fast` 启用，输入为 Feather V2，支持解析模式
`--permutations 0`。常规 `bgsmap run` 默认使用 `legacy` 引擎，支持原有格式和置换诊断。

## 在集群安装

在 wheel 所在目录执行。已有独立的 `bgsmap-lowmem 1.1.0` 时，先卸载它，再安装新版，
避免两个软件包同时管理 `bgsmap-lowmem` 命令。

```bash
conda activate gsmap
python -m pip uninstall -y bgsmap-lowmem
python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple \
  './bgsmap-1.2.0-py3-none-any.whl[all]'
bgsmap --version
bgsmap-lowmem --version
```

两条版本命令均应输出 `bgsmap 1.2.0`。也可在新版源码目录里用同一清华源安装
`'.[all]'`。本版通过 wheel 或源码安装，清华源提供依赖。

旧 HESTA 的 `prepare_hesta_bgsmap_lowmem.py` 会检查固定的 1.1.0 版本字符串。
本版使用下面的 Slurm 命令或直接执行 CLI，试跑输出指定新目录。

## 单表型、单切片 Slurm 试跑

把源码中的 `examples/hesta_lowmem.slurm` 复制到集群。替换四个输入位置和全新的结果目录，
在已激活安装环境的终端执行。

```bash
sbatch hesta_lowmem.slurm \
  /path/to/section_gene_marker_score.feather \
  /path/to/section.h5ad \
  /path/to/trait_burden.tsv.gz \
  /path/to/fast-scratch/gss-cache/section \
  /path/to/results-v1.2.0/trait/section

sacct -j JOB_ID --format=JobID,State,AllocCPUS,ReqMem,MaxRSS,Elapsed
```

示例沿用 8 CPU、100G 的试跑申请，安排 7 个 spot 工作线程和一个读取线程。
日志给出两遍扫描的各自耗时、其中读取块耗时及写表阶段。
任务生成四个结果表和 `run_metadata.json`。需绘图时可单独运行
`bgsmap plot --out /path/to/results-v1.2.0/trait/section --trait-name TRAIT`。

同一切片的所有表型传入同一个缓存目录。缓存以写完后原子重命名的方式发布，
重复调用 `cache-gss` 会校验并复用已有缓存。源路径、大小、修改时间或基因列不符时，
在新目录重建。缓存写入全部原始基因的 float32 数据，大小约为
`4 × Feather 基因行数 × Feather spot 列数` 字节。
若恰有 17,008 行和 563,390 列，需要约 35.7 GiB。旧日志中的 17,008 是匹配基因数，
因此实际大小应以 `cache-gss` 日志为准。首次转存需要读取 Feather 一遍，
每次分析仍需读取两遍连续缓存块。

## 无缓存直跑

磁盘空间不足时，省略缓存参数，直接读取 Feather。

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
bgsmap-lowmem run \
  --gss /path/to/section_gene_marker_score.feather \
  --h5ad /path/to/section.h5ad --annotation substructure \
  --burden /path/to/trait_burden.tsv.gz --gene-col gene_symbol --p-col Pvalue_Burden \
  --score-mode logp --duplicate-policy bonferroni --permutations 0 --seed 42 \
  --threads 7 --block-size 512 --projection-size 4096 --skip-plots \
  --out /path/to/results-v1.2.0/trait/section
```

- `--threads`（别名 `--n-threads`）控制同一表型内的 spot 块并发数。默认 1，
  也可用 `BGS_LOWMEM_THREADS` 设置。
- `--block-size` 控制每个计算块的 spot 数。默认 512，也接受
  `BGS_LOWMEM_SPOT_BLOCK`。试跑内存偏高时，可先减小块或减少线程。
- `--projection-size` 默认 4096，负责 Feather 投影与解压，对已建缓存无效。
  增大可降低宽表的重复元数据开销，同时增加瞬时内存。1.1.0 的 256 列块在
  563,390 spots 上每遍约 2,201 次投影，4,096 列对应约 138 次。
- `--skip-plots` 跳过绘图，也接受旧环境变量 `BGS_LOWMEM_SKIP_PLOTS=1`。
  H5AD 读取、结果表、排序和绘图有各自的额外开销，块和线程参数约束的是矩阵工作集。

## 统计一致性

1.2.0 在匹配的 spots 上做基因归一化，在匹配的基因及分层内计算与原版一致的解析
条件置换方差。区域包含区内平均和相对区外的竞争性结果。输出的 FWER 为 Holm，
FDR 为 BH。结果复用相同的解析 P 和校正函数。
流式 float64 归约和旧版 float32/不同顺序归约会产生舍入差异。
对于极近恒定且非零的基因，新版先中心化再累计二阶矩，以减少旧版原始平方和公式的消减误差。

小样本通过 LZ4 多 record batch Feather、H5AD 分类注释、不完全且重排的 spot
匹配、Bonferroni 去重、logp 和分层 binary、非有限/恒定基因 QC、
缓存复用和无效缓存拒绝测试。四张结果表均与旧引擎逐列比较。
使用 GitHub v1.0.0 源码和用户原始 1.1.0 wheel 也做了单独的同输入回归。

## 本地基准

在 Apple M2、Python 3.12、PyArrow 18.1、NumPy 1.26、局部 SSD 上运行，
两个数据集均是 81 个域、LZ4 单 record batch、float16 GSS、`--permutations 0`，
计入进程启动和结果表写入，跳过绘图。下表为每种配置一次端到端测量。
1.1.0 使用 256 spots/块，新版使用 512 spots/计算块及 4,096 列/Feather 投影。

| 模式 | 40,000 spots × 2,048 genes 时间 / 峰值 RSS | 16,384 spots × 17,008 genes 时间 / 峰值 RSS |
|---|---:|---:|
| 1.1.0 原始 wheel | 6.97 s / 154 MiB | 6.58 s / 266 MiB |
| 1.2.0 直读，1 个 worker | 1.94 s / 215 MiB | 2.62 s / 791 MiB |
| 1.2.0 直读，4 个 workers | 2.27 s / 238 MiB | 2.42 s / 1,003 MiB |
| 1.2.0 复用缓存，1 个 worker | 1.57 s / 154 MiB | 2.38 s / 347 MiB |
| 1.2.0 复用缓存，4 个 workers | 1.62 s / 210 MiB | 1.87 s / 808 MiB |
| 首次构建缓存 | 0.89 s / 186 MiB | 1.37 s / 627 MiB |

在 40,000-spot 样本上，单独测量 1.1.0 的两遍扫描分别为 3.42 s 和 3.49 s，
其中 Feather 读取/投影调用占 2.84 s 和 2.67 s。调用计时包含元数据处理、
块读取和解压，无法把这些因素分别归因。复用缓存后，在 17,008-gene 样本的
第二遍中，1 到 4 个工作线程将耗时从 0.91 s 降到 0.42 s。
采用较大的 Feather 投影和更多工作线程时，峰值 RSS 会增加。

复现时先在单独环境准备原始 `bgsmap 1.0.0` 和原 `bgsmap-lowmem 1.1.0` wheel，
再用安装新版的 Python 执行以下命令。输出包括各模式日志和 `metrics.json`。

```bash
python benchmarks/lowmem.py --out ./bench-wide \
  --baseline-python /path/to/bgsmap-v1.1-env/bin/python

python benchmarks/lowmem.py --out ./bench-gene-heavy --genes 17008 --spots 16384 \
  --baseline-python /path/to/bgsmap-v1.1-env/bin/python
```

本次未在 HESTA 的 563,390-spot 真数据与集群文件系统上运行，实际作业耗时、
内存峰值及队列资源可用性需由试跑确认。
