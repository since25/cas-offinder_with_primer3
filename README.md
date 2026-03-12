# hg38 脱靶分析与引物设计综合工具 (实验室V2版)

这是一个用于鉴定 CRISPR sgRNA 潜在脱靶位点并自动为其设计验证引物的综合流程工具。
本工具底层使用 [Cas-OFFinder](https://github.com/snugel/cas-offinder) 在全基因组范围内进行快速脱靶检索，使用 `primer3-py` 提取侧翼序列并自动设计引物，结合 `pyranges` 提供 GTF 基因注释，最后通过 `streamlit` 提供直观的可视化 Web 界面。

## 🌟 核心特性

- **多核并行与 GPU 加速**：支持多线程查表，如果有 OpenCL 环境，Cas-OFFinder 可自动开启 GPU 计算。
- **Cas-OFFinder V3 支持**：除常规错配（Mismatches）外，支持最新的 **DNA 凸起 (DNA Bulge)** 和 **RNA 凸起 (RNA Bulge)** 容错检索。
- **靶向引物智能设计**：针对每个脱靶位点自动提取上游和下游侧翼序列，规避靶点区域 (额外 20bp 缓冲)，动态设计最优的扩增引物，并支持失败后的降级策略（放宽 Tm/GC 条件）。
- **Docker 极速部署**：提供高度精简的 Docker 环境，只需克隆代码即可拥有所有底层编译好的环境。
- **强大的 UI 界面**：通过 Streamlit 实现，支持调节靶点序列、PAM、错配度、扩增产物长度范围 (Amplicon Size) 等参数，并提供批量处理 (CSV) 能力。

## 📁 目录结构

```text
offtarget_primer_tool/
├── README.md              # 本说明文档
├── Dockerfile             # 极简 Docker 构建脚本
├── pyproject.toml         # Python 依赖及打包配置
├── setup.py               # 安装入口
├── build_linux/           # 存放 Linux 预编译好的 cas-offinder 执行文件
├── app/
│   └── streamlit_app.py   # Streamlit 网页端入口
├── data/
│   └── hg38/
│       ├── hg38.fa            # (需自行准备) 人类参考基因组 FASTA
│       └── annotation.gtf.gz  # (需自行准备) 基因注释文件
├── runs/
│   └── .cache/            # (自动生成) 历史查询缓存数据
└── src/
    └── otp/               # 核心层代码库 (脱靶检索、引物分析、并行管线、可视化等)
```

## 🛠 数据准备

为了使系统正常运行，您需要准备 hg38 基因组数据：
1. 下载 `hg38.fa` 并放入 `data/hg38/` 目录。
2. 生成或准备 `hg38.fa.fai` 索引文件。
3. (可选) 下载 GTF 基因注释文件并放入 `data/hg38/`，用于后续注释。

推荐目录结构：

```text
data/
└── hg38/
    ├── hg38.fa
    ├── hg38.fa.fai
    └── annotation.sorted.gtf.gz
```

### 下载 hg38 FASTA

可以从 UCSC 下载整套 hg38 主染色体 FASTA：

```bash
mkdir -p data/hg38
cd data/hg38
curl -L -o hg38.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

如果您习惯使用 `wget`：

```bash
mkdir -p data/hg38
cd data/hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

### 生成 FASTA 索引 `hg38.fa.fai`

推荐使用 `samtools faidx`：

```bash
cd data/hg38
samtools faidx hg38.fa
```

如果系统里没有 `samtools`，macOS 可以先安装：

```bash
brew install samtools
```

本项目底层使用 `pyfaidx` 读取 FASTA；如果缺少 `.fai`，多数情况下也会在首次读取时自动生成，但建议提前手工生成，便于排查问题。

### 下载基因注释 GTF

注释文件不是必须的；缺失时流程仍可运行，只是结果中的基因注释字段会被标记为跳过。

推荐直接下载 GENCODE 的 hg38 注释，然后按项目默认文件名保存为 `annotation.sorted.gtf.gz`：

```bash
cd data/hg38
curl -L -o annotation.sorted.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
```

如果您下载的是别的文件名，也可以手工重命名：

```bash
mv gencode.v47.annotation.gtf.gz annotation.sorted.gtf.gz
```

### 快速检查数据是否就绪

```bash
ls -lh data/hg38/hg38.fa data/hg38/hg38.fa.fai
```

如果还准备了注释文件：

```bash
ls -lh data/hg38/annotation.sorted.gtf.gz
```

> **注意**：由于基因组文件体积巨大（~3GB），代码仓库内的 `.gitignore` 已配置忽略 `data/` 目录。在部署到新服务器时，请务必手动上传这些大文件。

## 🚀 安装与运行 (极简 Docker 部署)

强烈推荐使用 Docker 运行，彻底告别环境和依赖报错。在克隆本仓库并放好数据文件后：

### 1. 构建镜像
由于仓库中已经包含了编译好的 Linux 格式 Cas-OFFinder 二进制文件 (`build_linux/cas-offinder`)，Docker 构建速度极快：
```bash
docker build -t cas_primer_app .
```

### 2. 运行容器页面
使用如下命令，启动容器并将您的本地基因组数据挂载进系统：
```bash
docker run -p 8501:8501 -v "$(pwd)/../data:/app/data" --name cas_primer_app_container cas_primer_app
```
启动后，直接在浏览器中访问 `http://localhost:8501` 即可使用网页界面。


## 💻 开发与本地部署

如果您希望在 macOS 或本地环境直接开发：

```bash
# 激活虚拟环境
python3 -m venv .venv
source .venv/bin/activate

# 安装本体依赖
pip install -e .

# 启动网页端
streamlit run app/streamlit_app.py
```

## 🧭 常见使用入口

### 1. Streamlit 网页端
适合手工输入单条 sgRNA，或上传批量参数表后直接点击运行：

```bash
source .venv/bin/activate
streamlit run app/streamlit_app.py
```

默认参数说明：
- `Flank Length` 默认 `500`，表示脱靶位点上下游各截取 `500 bp`
- Primer3 在目标区域两侧还会额外避让 `20 bp`
- Streamlit 内部会自动复用当前 Python 解释器调用 `otp.pipeline`

### 2. 命令行批量模式：从 sgRNA 参数表开始跑全流程
适合输入是“待搜索参数表”，也就是还没有 Cas-OFFinder 结果，需要从 sgRNA/PAM/mismatch 开始做全流程搜索和引物设计。

支持 `CSV` 和 `Excel`：

```bash
source .venv/bin/activate
python -m otp.pipeline --batch input.xlsx --out runs/batch_run --threads 8
```

批量输入表至少应包含：
- `spacer`
- 可选：`pam`, `mismatches`, `dna_bulge`, `rna_bulge`, `flank`, `amplicon_min`, `amplicon_max`

单条命令行模式：

```bash
source .venv/bin/activate
python -m otp.pipeline \
  --spacer GAGTCCGAGCAGAAGAAGA \
  --pam NGG \
  --mismatches 3 \
  --out runs/single_query
```

### 3. 命令行重算模式：已有脱靶结果 Excel，直接重设计引物
适合输入已经是 Cas-OFFinder 或旧脚本导出的脱靶位点表，不需要重新跑 Cas-OFFinder，只想基于现有位点重新设计引物。

正式入口：

```bash
source .venv/bin/activate
python -m otp.redesign \
  --input your_offtargets.xlsx \
  --out runs/your_offtargets_redesign
```

安装为脚本后也可以直接用：

```bash
source .venv/bin/activate
otp-redesign \
  --input your_offtargets.xlsx \
  --out runs/your_offtargets_redesign
```

这类输入表至少应包含以下字段：
- `crRNA`
- `DNA`
- `Chromosome`
- `Position`
- `Direction`
- `Mismatches`
- `Bulge Size`

可选字段：
- `#Bulge type` 或 `Bulge type`
- `eff_pred`

示例：

```bash
source .venv/bin/activate
python -m otp.redesign \
  --input sgRNA-ANGPTL3-EXON1-CAA4.xlsx \
  --out runs/sgRNA-ANGPTL3-EXON1-CAA4_legacy_export
```

### 4. 什么时候走哪个入口
- 如果输入是 sgRNA、PAM、mismatch 这些查询参数，用 `otp.pipeline --batch`
- 如果输入已经是脱靶位点 Excel，用 `otp.redesign --input`
- 如果想在浏览器里操作，用 Streamlit

## 📊 输出结果说明

流程结束后，会在 `runs/` 目录下生成 `results.csv` 与 `results.xlsx`。

`results.xlsx` 默认包含以下 sheet：
- `offtargets`：脱靶位点基础信息
- `primers`：引物设计结果与坐标
- `legacy_primers`：兼容旧版实验室导出格式的中文字段表
- `summary`：结果汇总

主要字段说明：
- `chrom`, `pos0`, `strand`: 脱靶位点的染色体与坐标。
- `mismatches`, `bulge_type`, `bulge_size`: 编辑距离及凸起判定。
- `primer_left_seq`, `primer_right_seq`, `amplicon_size`: 软件挑选出的最优上、下游引物序列及预估的产物长度。
- `covers_offtarget`: 验证标记（为 True 说明扩增产物准确包含了发生错配的脱靶核心区）。

坐标相关字段补充：
- `primer_left_pos_in_flank0` / `primer_right_pos_in_flank0`：引物在内部截取的 flank 序列中的 `0-based` 位置
- `flank_start0` / `flank_end0`：内部 flank 序列对应的基因组区间，默认是脱靶位点上下游各 `500 bp`
- `primer_left_genome0` / `primer_right_genome0`：引物在基因组上的 `0-based` 坐标

`legacy_primers` 中的旧格式字段说明：
- `正向引物结合位点`、`反向引物结合位点`：相对于 `位点上下游1000bp序列` 的位置
- `位点上下游1000bp序列`：按脱靶位点上下游各 `1000 bp` 提取
- `扩增子序列`：按当前左右引物基因组坐标提取的扩增片段

## 💡 缓存机制
系统会自动将输入参数（Spacer、mismatch 等参数）的 Hash 值记录到 `runs/.cache/`，当下一次参数完全相同时，直接秒速恢复历史运算结果，省去 Cas-OFFinder 漫长的搜索时间。
