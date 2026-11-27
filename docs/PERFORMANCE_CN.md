# bgsmap 1.3.0 大数据性能与内存指南

bgsmap 在匹配的基因和空间 spots 上计算 burden 富集。计算与调度器相互独立，
可在个人电脑、工作站及集群计算节点通过同一条命令使用。
输入 GSS 为 gene × spot，支持 Feather V2、Parquet、TSV 和 CSV，
空间坐标与域标签来自任意匹配的 H5AD 或坐标表。

## 安装与普通入口

在下载的 wheel 所在目录执行。旧的独立低内存扩展与新版共用命令名，
安装新版前先卸载该扩展。

```bash
python -m pip uninstall -y bgsmap-lowmem
python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple \
  './bgsmap-1.3.0-py3-none-any.whl[all]'
bgsmap --version
```

源码安装使用 `python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple '.[all]'`。
分发包名和导入名均为 `bgsmap`。保留的 `bgsmap-lowmem` 命令默认显式选择 fast 引擎。

```bash
bgsmap run \
  --gss /path/to/GSS.feather \
  --h5ad /path/to/spatial.h5ad --annotation domain \
  --burden /path/to/burden.tsv.gz --gene-col gene --p-col p \
  --out /path/to/results
```

默认 `--engine auto` 在解析模式 `--permutations 0` 下选择分块引擎，
自动选择线程及矩阵块。选择 `--engine legacy` 可运行原内存引擎做对照。
请求置换诊断或读取旧 Feather V1 文件时，auto 选择原内存引擎。fast 目前支持解析模式。
`quick` 在上游 gsMap 生成 GSS 后走同一入口，其训练阶段的资源由 gsMap 管理。

改变文件格式、注释列、基因/P 列、组织类型和 spot 命名方式均无需修改程序。
若使用坐标表，以 `--coords spatial_coords.tsv` 替换 H5AD 和注释参数，
坐标表包含 `spot,x,y,domain` 列。仍按完整 spot ID 匹配并保持 GSS 顺序。
基因映射、重复测试校正和匹配协变量沿用既有选项。

## 按可用资源运行

```bash
bgsmap run \
  --gss /path/to/GSS.parquet --coords /path/to/spatial_coords.tsv \
  --burden /path/to/burden.tsv.gz \
  --threads 8 --memory-limit 4G --skip-plots \
  --out /path/to/results
```

- `--threads auto` 默认使用最多四个 spot 工作线程，Linux 上参考进程 CPU 亲和性。
  指定整数可使用更多核，每个工作线程处理同一表型内的不同 spot 块。
- `--memory-limit 4G` 是分析工作集估算预算。默认取系统/容器可用内存一半，
  上限 4 GiB。估算包括矩阵、域聚合及结果表预留空间。
- `--block-size` 指定每个计算块的 spots，自动候选值为 512。
  `--projection-size` 指定每次 Feather/Parquet 读取的列数，自动候选值为 4096。
  自动设置会按预算缩小。显式值保留，在估算超预算时给出错误。
- `--skip-plots` 可把分析和绘图分开运行。之后执行
  `bgsmap plot --out /path/to/results` 生成 HTML 和矢量 PDF。

内存预算用于规划分配。实际 RSS 还受到 Python/Arrow 分配器、文件映射和绘图影响。
文本导入还有解析器内存开销。需要严格内存上限时，由运行环境限制进程内存。
日志打印采用的资源计划，`run_metadata.json` 的 `low_memory` 字段保存相同计划。

旧环境变量 `BGS_LOWMEM_THREADS`、`BGS_LOWMEM_SPOT_BLOCK` 和
`BGS_LOWMEM_SKIP_PLOTS` 继续有效。相应命令行参数优先。

## 一份 GSS 的重复分析

```bash
bgsmap run \
  --gss /path/to/GSS.feather --h5ad /path/to/spatial.h5ad --annotation domain \
  --burden /path/to/trait_A.tsv.gz \
  --gss-cache /path/to/fast-disk/sample-cache \
  --memory-limit 4G --threads 8 --out /path/to/results/trait_A

bgsmap run \
  --gss /path/to/GSS.feather --h5ad /path/to/spatial.h5ad --annotation domain \
  --burden /path/to/trait_B.tsv.gz \
  --gss-cache /path/to/fast-disk/sample-cache \
  --memory-limit 4G --threads 8 --out /path/to/results/trait_B
```

缓存目录不存在时自动构建，之后复用。也可提前执行
`bgsmap cache-gss --gss GSS.feather --out sample-cache --memory-limit 4G`。
缓存只保存 float32 原始 GSS，分析时根据各自匹配的 spots 和基因进行 QC 与归一化，
因此可复用于不同表型、注释和匹配基因集合。既有 1.2.0 缓存可继续校验和读取。
缓存大小约为 `4 × 原始基因数 × 原始 spot 数` 字节。
缓存保存源文件路径、大小、修改时间和基因列，源文件变化或缓存被截断会拒绝复用。
完成后的目录通过原子重命名发布，失败时清理临时目录。

Feather/Parquet 采用 spot-major 缓存，计算阶段按连续块读取。
TSV/CSV 按基因行分块解析，采用 gene-major 缓存，并在分析时只复制选定矩阵块。
文本没有指定持久缓存时，程序在系统临时目录建立缓存，并在成功或失败后清理。
可用 `--temp-dir /path/to/large-scratch` 选择有足够磁盘空间的位置。
首次文本转存需要额外时间，重复分析可复用 `--gss-cache` 省去转存。

## 定位耗时

`run_metadata.json` 的 `performance` 记录源数据准备、对齐、基因 QC、
spot 计算、写表、绘图和总耗时。两遍扫描还记录读取耗时、工作线程累计计算耗时、
块数量及线程数。读取与计算会重叠，累计工作线程时间可能超过扫描的墙钟时间。

若读取占比较高，可在快速本地磁盘复用缓存或调整投影大小。
若 spot 计算占比较高，可比较 1、2、4、8 线程及不同块大小的实测结果。
若写表或绘图占比较高，可查看对应阶段时间并单独运行绘图。

## 统计与验证范围

分块引擎使用与原版相同的 burden 权重、分层内条件置换方差、
单侧正态尾部解析 P、Holm FWER 和 BH FDR。域结果包含区内均值与相对区外的
竞争性统计量。基因归一化在当前分析匹配的全部 spots 上完成。
float64 流式归约与原内存引擎的 float32 归约存在舍入差异。

测试覆盖直接 CLI 和 Python 入口、不同输入格式与压缩方式、多 record batch /
多 row group、不同线程与块、重排及不完全匹配的 spots、单域与多域、
logp/binary 权重、协变量分层、Bonferroni 去重、恒定及非有限基因、
预算调整、缓存复用、无效缓存拒绝、临时缓存异常清理及绘图。
原版解析统计的假设仍适用，性能优化本身无法验证任意真数据的零假设校准。

## 跨格式合成基准

以下为 Apple M2、Python 3.12、NumPy 1.26、PyArrow 18.1、局部 SSD 上的单次测量，
各配置依次运行。时间包含进程启动与写表，跳过绘图。
auto 使用默认最多四线程和 1 GiB 估算预算。缓存列为已有缓存后的运行，
首次缓存构建耗时单列。缓存大小、计算块和元数据均随实际矩阵变化。

| 输入与规模 | 对照耗时 / 峰值 RSS | auto 耗时 / 峰值 RSS | 4 线程复用缓存耗时 / 峰值 RSS | 首次建缓存 |
|---|---:|---:|---:|---:|
| Zstd Feather，50,000 spots × 768 genes，18 域 | 1.1.0 为 9.21 s / 154 MiB | 2.35 s / 169 MiB | 1.41 s / 150 MiB | 0.91 s |
| Parquet，4,096 spots × 2,048 genes，12 域 | 原内存引擎为 0.89 s / 420 MiB | 0.84 s / 328 MiB | 0.54 s / 179 MiB | 0.58 s |
| 压缩 TSV，4,096 spots × 1,024 genes，6 域 | 原内存引擎为 0.73 s / 296 MiB | 1.23 s / 298 MiB | 0.52 s / 145 MiB | 1.05 s |

Zstd Feather 对照使用用户提供的原始 1.1.0 wheel 及 GitHub v1.0.0 基础包。
其他格式的对照使用 1.3.0 中的 legacy 引擎，以便统一跳过绘图。
四张输出表逐列通过数值回归比较。TSV 首次临时转存的开销在这个规模下使单次运行变慢。
实际加速取决于格式、存储和矩阵形状，增加线程也可能受读取或内存带宽限制。

复现时为第一个命令提供安装原始 1.1.0 wheel 与 bgsmap 1.0.0 的 Python。
运行命令的主 Python 安装当前版本。

```bash
python benchmarks/lowmem.py --out bench-wide \
  --genes 768 --spots 50000 --domains 18 --compression zstd \
  --baseline-python /path/to/bgsmap-1.1-env/bin/python

python benchmarks/lowmem.py --out bench-parquet \
  --genes 2048 --spots 4096 --domains 12 --format parquet --baseline-engine legacy

python benchmarks/lowmem.py --out bench-text \
  --genes 1024 --spots 4096 --domains 6 --format tsv.gz --baseline-engine legacy
```

基准输出逐模式的日志、结果表和 `metrics.json`。源数据生成在计时之外。
完整 563,390-spot HESTA 以及其他平台的超大真数据、网络文件系统仍需实测。
