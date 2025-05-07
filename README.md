# amplysis
![version](https://img.shields.io/badge/version-1.3.2-blue.svg)
![license](https://img.shields.io/badge/license-GPL--3.0-green.svg)
![platform](https://img.shields.io/badge/platform-R-75AADB.svg)
![status](https://img.shields.io/badge/status-active-brightgreen.svg)
![CRAN](https://img.shields.io/badge/CRAN-Submitting-orange.svg)

**An R package for rapid analysis of 16S amplicon sequencing data**
> **16S 扩增子测序数据分析 R 包**

An R package for 16S rRNA gene amplicon sequencing data that integrates data preprocessing, analysis, and visualization methods. These methods include microbial composition analysis, α-diversity analysis, β-diversity analysis, differential analysis, correlation analysis, and network analysis, among others.

> 一个用于 16S rRNA 基因扩增子测序数据，集成了数据预处理、分析与可视化方法的 R 包。这些方法包括微生物组成分析、α 多样性分析、β 多样性分析、差异分析、关联分析和网络分析等。

<br>

## Background 背景

The rapid development of bioinformatics tools has enabled researchers to perform complex 16S rRNA gene amplicon sequencing analysis. However, for those without a bioinformatics background, existing tools and R packages can be complex and difficult to use. The amplysis project was created to address this issue, providing an accessible solution for researchers to easily conduct data analysis. 

> 生物信息学工具的快速发展，使得研究者能够进行复杂的 16S rRNA 基因扩增子测序数据分析。然而，对于缺乏生物信息学背景的研究者来说，现有的工具和 R 包操作复杂，使用门槛较高。amplysis 集合 R 包项目应运而生，旨在为研究者提供一个易于使用的数据分析解决方案。  

<br>

## Install 安装
### 1. Install the latest version from GitHub (recommended)
> 1. 从 Github 上安装最新版本（推荐）

安装 R 包 `devtools`：

```
install.packages("devtools")
library(devtools)
```
通过 `devtools` 安装 R 包 `amplysis`：
```
devtools::install_github("min-perilla/amplysis")
```

### 2. Install from a local `.tar.gz` source package
> 2. 本地 `.tar.gz` 源码包安装

Go to the [Releases page](https://github.com/min-perilla/amplysis/releases) to download the latest `amplysis_*.tar.gz` release package and note its save location. Then use the following code in R to install the local source package (replace the path with your actual file location)
> 前往 [Releases 页面](https://github.com/min-perilla/amplysis/releases) 下载最新发布的 `amplysis_*.tar.gz` 安装包，并记住其保存路径。使用以下代码在 R 中安装本地源码包（请将路径替换为你实际的文件位置）

Windows users can install by entering the following code in RStudio:
> Windows 用户在 RStudio 中输入下列代码进行安装:
```
install.packages("C:/path/to/your/amplysis_X.X.X.tar.gz", repos = NULL, type = "source")
```

macOS users can install by entering the following code in RStudio:
> macOS 用户在 RStudio 中输入下列代码进行安装:
```
install.packages("/Users/path/to/your/amplysis_X.X.X.tar.gz", repos = NULL, type = "source")
```

### 3. Install from CRAN
> 3. 从 CRAN 安装

> [!IMPORTANT]
>
> Submitting to CRAN, stay tuned...
>
> 正在向 CRAN 提交中，敬请期待...

```
install.packages("amplysis")
```

<br>

## Usage 使用方法

Click [here](https://github.com/min-perilla/amplysis/releases/download/latest/example_data.zip) to download the example data. 
> 点击[这里](https://github.com/min-perilla/amplysis/releases/download/latest/example_data.zip) 下载示例数据文件。

After clicking the download link, you will get a compressed file named `example_data.zip`. Extract the zip file and open the folder, which contains a test script named `example.R` and a `data` directory (containing six types of example data files). 

You can open the `example.R` file in RStudio, select all the code, and click `Run` to execute the test. Alternatively, you may copy the sample code below into `RStudio` for testing:

> 点击下载链接后，会得到一个名为 `example_data.zip` 的压缩包。解压并打开文件夹，里面包含了一个名为 `example.R` 的测试脚本和一个 `data` 文件夹（内含 6 种类型的示例数据文件）。
>
> 你可以在 `RStudio` 中打开 `example.R` 文件，全选代码并点击 `Run` 运行测试。或者，您也可以将下方的示例代码复制到 `RStudio` 中进行测试：

```
# example.R
# Clear all variables | 清除所有变量
rm(list = ls())

# Load amplysis package | 加载 amplysis
library(amplysis)

# Set working directory | 设置工作目录
set_wd()


# ------------------------------------------------------------------------------
# Load example data | 加载示例数据
otu = read_data("data/otu.csv")            # OTU table | 特征表
tax = read_data("data/tax.csv")            # Taxonomy table | 分类表
metadata = read_data("data/metadata.csv")  # Sample metadata | 样本元数据
rep = read_data("data/rep_seqs.csv")       # Representative sequences | 代表性序列
env = read_data("data/env.csv")            # Environmental factors | 环境因子
tree = read_data("data/tree_rooted.nwk")   # Phylogenetic tree | 系统发育树


# ------------------------------------------------------------------------------
# Data preprocessing: Taxonomy table | 数据预处理：分类表
# Split taxonomy table into columns | 分类表数据分列
tax = tax_separate(tax = tax, index = 2, delim = "; ")      
# Remove prefixes from taxonomy table | 分类表去前缀
tax = tax_trim_prefix(tax = tax, index = c(2:8), length = 3)
# Repair taxonomy table names | 分类表物种注释信息优化
tax = tax_names_repair(tax = tax, column_to_check = 7, column_to_add = 3)

# Data preprocessing: Rarefaction | 数据预处理：数据抽平
otu_tax = data_rarefy(otu, method = "phyloseq", tax_table = tax)

# Separate OTU table and taxonomy table | 拆分特征表和分类表
otu = otu_tax[["otu"]]
tax = otu_tax[["tax"]]

# Align representative sequences file | 对齐代表性序列文件
otu_rep = merge(x = otu, y = rep, by = "#OTU ID", all.x = T, sort = F)

# Separate OTU table and representative sequences | 拆分特征表和代表性序列文件
rep = otu_rep[, c(1, ncol(otu_rep))]
otu = otu_rep[, -c(ncol(otu_rep))]


# ------------------------------------------------------------------------------
# Data Analysis & Visualization | 数据分析与可视化

# Stacked bar plot analysis (Phylum level) | 物种堆叠图分析（门水平）
data_sta_p = stackbar(otu, tax, metadata, tax_cla = "phylum", group1 = "group", group2 = "group2", row_n = 8)
# Visualization | 可视化
stackbar_plot(data_sta_p, tax_cla = "phylum", title_legend = "Top 8 Phyla")


# Stacked bar plot analysis (Genus level) | 物种堆叠图分析（属水平）
data_sta_g = stackbar(otu, tax, metadata, tax_cla = "genus", group1 = "group", group2 = "group2", row_n = 20)
# Visualization | 可视化
stackbar_plot(data_sta_g, tax_cla = "genus", title_legend = "Top 20 Genera")


# Chord diagram (Phylum level) | 弦图（门水平）
data_chord = chord(otu, metadata, tax, tax_cla = "phylum", group = "group2", row_n = 8)
# Visualization | 可视化
chord_plot(data_chord)


# Venn diagram | 韦恩图
data_venn = venn(otu, metadata, group = "group2")
# Visualization | 可视化
venn_plot(data_venn)


# Upset plot | 集合图
data_upset = Upset(otu, metadata, group = "group2")
# Visualization | 可视化
Upset_plot(data_upset)


# Boxplot (Alpha diversity analysis) | 箱线图（Alpha 多样性分析）
data_alpha = alpha(otu, metadata, group = "group2", tree = tree)
# Visualization | 可视化
alpha_plot(data_alpha)


# Principal Component Analysis (PCA) | 主成分分析（PCA）
data_pca = pca(otu, metadata, group = "group2")
# Visualization | 可视化
pca_plot(data_pca)


# Principal Coordinates Analysis (PCoA) | 主坐标分析（PCoA）
data_pcoa = pcoa(otu, metadata, group = "group2")
# Visualization | 可视化
pcoa_plot(data_pcoa)


# Non-metric Multidimensional Scaling (NMDS) | 非度量多维尺度分析（NMDS）
data_nmds = nmds(otu, metadata, group = "group2")
# Visualization | 可视化
nmds_plot(data_nmds)


# Redundancy Analysis (RDA) | 冗余分析（RDA）
data_rda = RDA(otu, env, metadata, group = "group2")
# Visualization | 可视化
RDA_plot(data_rda)


# Canonical Correspondence Analysis (CCA) | 典型对应分析（CCA）
data_cca = CCA(otu, env, metadata, group = "group2")
# Visualization | 可视化
CCA_plot(data_cca)


# Heatmap | 热图
data_heatmap = heatmap(otu, tax, metadata, tax_cla = "genus", group1 = "group", group2 = "group2", row_n = 30)
# Visualization | 可视化
heatmap_plot(data_heatmap, fontsize_col = 14, file_height = 10, file_width = 12)


# Co-occurrence network analysis | 共现性网络分析
# Visualization using Gephi is recommended | 推荐使用 Gephi 软件可视化
data_net = network(otu, tax, metadata, tax_cla = "genus")
data_net

```

Partial example figures:
> 部分示例图：

![image](https://github.com/user-attachments/assets/1a2e6ba6-6e1f-4846-97f7-81332c1d14b9)

<br>

## Additional Information 补充说明
### 1. Preparation of Data Files
> 数据文件的准备

Please refer to the example files in the `test data/example data` directory on GitHub to prepare your data accordingly for analysis. For more detailed information about the data files, please refer to this article (DOI will be provided once available).
> 请参考目录 `/test data/example data` 下的示例文件，以相应格式准备好数据文件进行分析。关于数据文件的更多详细说明，可以参考此文章（DOI 可用后将提供）。
<br>

### 2. Using `read_data()` to Read `.biom` Files
> 使用函数 `data_rarefy()` 读取 .biom 文件

The `read_data()` function supports reading `JSON` format `.biom` files. If you already have a `.biom` file (such as one exported from `QIIME 2`, which may be in `HDF5` format), you can convert it to `JSON` format using the following methods.
> `read_data()` 函数支持读取 `JSON` 格式的 `.biom` 文件。如果你已经有一个 `.biom` 文件（如 `QIIME 2` 导出的 `.biom` 文件，可能是 `HDF5` 格式），你可以通过以下方法将其转换为 `JSON` 格式。

#### Windows：

1. Check if `biom` is already installed:
> 1. 检查 `biom` 是否已经安装：
   ```bash
   biom --version
   ```

   If not installed, use the following command to install the `biom-format` tool:
   > 如果未安装，可以使用以下命令安装 `biom-format` 工具：

   ```bash
   pip install biom-format
   ```

2. Convert the `HDF5` format `.biom` file to `JSON` format:
> 2. 使用以下命令将 `HDF5` 格式的 `.biom` 文件转换为 `JSON` 格式：

   ```bash
   biom convert -i input_table.biom -o output_table_json.biom --to-json --table-type "OTU table"
   ```
<br>

#### macOS：

1. Check if `biom` is already installed:
> 1. 检查 `biom` 是否已经安装：

   ```bash
   biom --version
   ```

   If not installed, use the following command to install the `biom-format` tool:
   > 如果未安装，可以使用以下命令安装 `biom-format` 工具：

   ```bash
   pip install biom-format
   ```

2. Convert the `HDF5` format `.biom` file to `JSON` format:
> 2. 使用以下命令将 `HDF5` 格式的 `.biom` 文件转换为 JSON 格式：

   ```bash
   biom convert -i input_table.biom -o output_table_json.biom --to-json --table-type "OTU table"
   ```
<br>

#### Warning When Reading `.biom` Files in `HDF5` Format
> 关于读取 `HDF5` 格式的 `.biom` 文件的警告

If you attempt to read an `HDF5`-format `.biom` file (such as one exported from `QIIME 2`), the `read_data()` function may produce the following warning:
> 如果你尝试读取 `HDF5` 格式的 `.biom` 文件（如 `QIIME 2` 导出的 `.biom` 文件），`read_data()` 函数可能会出现以下警告：

```R
Warning messages:
1: In strsplit(conditionMessage(e), "\n") :
  unable to translate 'lexical error: invalid char in json text.
                                       <89>HDF                        (right here) ------^
' to a wide string
2: In strsplit(conditionMessage(e), "\n") : input string 1 is invalid
```
Although these warnings may appear, the file can still be read correctly. Rest assured, this issue will be resolved in a future update.
> 虽然会出现这些警告，但文件是可以正常读取的。请放心，未来我们会解决此问题。
<br>

### 3. About the data_rarefy() function 
> 关于函数 data_rarefy() 的说明

The R package provides a function `data_rarefy()` for data rarefaction. It integrates the `rrarefy()` function from the `vegan` package or the `rarefy_even_depth()` function from the `phyloseq` package for data rarefaction.  
  
Note that for the `rarefy_even_depth()` function in the `phyloseq` package, both replacement sampling and non-replacement sampling methods can be chosen. In replacement sampling, each time an element is selected from the sample pool, it is returned, keeping the sample pool unchanged, allowing some elements to be selected multiple times. This method can speed up computation and reduce memory usage, especially when handling large datasets. However, it may result in some OTU or ASV counts exceeding their original values, which could affect the accuracy of the analysis. 
  
In contrast, non-replacement sampling removes the selected sample after each draw, ensuring that each sample is selected only once, thereby ensuring OTU/ASV counts do not exceed their original values. The `rarefy_even_depth()` function from the `phyloseq` package enables replacement sampling by default (parameter `replace = TRUE`), which may have a slight impact on the resulting analysis. To ensure accuracy, the `data_rarefy()` function from the `amplysis` package uses non-replacement sampling for data rarefaction (parameter `replace = FALSE`).  

> R 包提供了一个函数 `data_rarefy()` 用于数据稀释（rarefaction）。它集成了 `vegan` 包中的 `rrarefy()` 函数或 `phyloseq` 包中的 `rarefy_even_depth()` 函数进行数据稀释。  
> 
> 需要注意的是，对于 `phyloseq` 包中的 `rarefy_even_depth()` 函数，可以选择替代抽样（replacement sampling）和非替代抽样（non-replacement sampling）方法。在替代抽样中，每次从样本池中选取一个元素后，该元素会被放回，使样本池保持不变，因此某些元素可能会被多次选中。这种方法可以加快计算速度并减少内存使用，特别适用于处理大型数据集。但它可能导致某些 OTU 或 ASV 的计数超过原始值，从而影响分析结果的准确性。
>   
> 相比之下，非替代抽样在每次抽样后移除已选中的样本，确保每个样本仅被选中一次，从而保证 OTU/ASV 的计数不会超过其原始值。`phyloseq` 包中的 `rarefy_even_depth()` 函数默认启用替代抽样（参数 `replace = TRUE`），这可能对最终分析结果产生轻微影响。为确保结果准确，`amplysis` 包中的 `data_rarefy()` 函数采用非替代抽样进行数据稀释（参数 `replace = FALSE`）。

<br>

## Maintainers 项目主要负责人
[@min-perilla](https://github.com/min-perilla)

<br>

## Contributing 贡献
[@min-perilla](https://github.com/min-perilla)

<br>

## License 许可证
The project is licensed under **GPL (>= 3)**. For more details, please refer to the **[LICENSE](https://github.com/min-perilla/amplysis/blob/main/LICENSE)** file.
> 该项目采用 **GPL (>= 3)** 许可证，详情请参阅 **[LICENSE](https://github.com/min-perilla/amplysis/blob/main/LICENSE)** 文件。
