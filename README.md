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
1. 下载 `hg38.fa` 以及它的索引文件 `hg38.fa.fai`，放入 `data/hg38/` 目录下。
2. (可选) 下载 GTF 基因注释文件 (如 `annotation.gtf.gz`) 放入 `data/hg38/`，用于后续 PyRanges 基因位点注释。

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
python3 -m venv venv
source venv/bin/activate

# 安装本体依赖
pip install -e .

# 启动网页端
streamlit run app/streamlit_app.py
```

### 命令行批量模式
```bash
python -m otp.pipeline --batch input.csv --out runs/batch_run/ --threads 8
```

## 📊 输出结果说明

流程结束后，会在 `runs/` 目录下生成清晰的 `results.csv` 与 `results.xlsx`（按标签页排版），主要包含：
- `chrom`, `pos0`, `strand`: 脱靶位点的染色体与坐标。
- `mismatches`, `bulge_type`, `bulge_size`: 编辑距离及凸起判定。
- `primer_left_seq`, `primer_right_seq`, `amplicon_size`: 软件挑选出的最优上、下游引物序列及预估的产物长度。
- `covers_offtarget`: 验证标记（为 True 说明扩增产物准确包含了发生错配的脱靶核心区）。

## 💡 缓存机制
系统会自动将输入参数（Spacer、mismatch 等参数）的 Hash 值记录到 `runs/.cache/`，当下一次参数完全相同时，直接秒速恢复历史运算结果，省去 Cas-OFFinder 漫长的搜索时间。
