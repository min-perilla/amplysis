# amplysis
![version](https://img.shields.io/badge/version-1.3.0-blue.svg)
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
### Install from CRAN
> 从 CRAN 安装

> [!IMPORTANT]
>
> Submitting to CRAN, stay tuned...
>
> 正在向 CRAN 提交中，敬请期待...

```
install.packages("amplysis")
```

### Install the latest version from GitHub (recommended)

> 从 Github 上安装最新版本（推荐）

安装 R 包 `devtools`：

```
install.packages("devtools")
library(devtools)
```
通过 `devtools` 安装 R 包 `amplysis`：
```
devtools::install_github("min-perilla/amplysis")
```

<br>

## Usage 使用方法
[示例数据下载链接](https://github.com/min-perilla/amplysis/releases/download/Latest/example.data.zip)
```
# 清除所有变量 | Clear all variables
rm(list = ls())

# 加载 amplysis | Load amplysis package
library(amplysis)

# 设置工作目录 | Set working directory
set_wd()

# 加载示例数据 | Load example data
otu = read_data("otu.csv")            # 特征表 | OTU table
tax = read_data("tax.csv")            # 分类表 | Taxonomy table
metadata = read_data("metadata.csv")  # 样本元数据 | Sample metadata
rep = read_data("rep_seqs.csv")       # 代表性序列 | Representative sequences
env = read_data("env.csv")            # 环境因子 | Environmental factors
tree = read_data("tree_rooted.nwk")   # 系统发育树 | Phylogenetic tree

# 数据预处理：分类表 | Data preprocessing: Taxonomy table
tax = tax_separate(tax = tax, index = 2, delim = "; ")                     # 分类表数据分列 | Split taxonomy table into columns
tax = tax_trim_prefix(tax = tax, index = c(2:8), length = 3)               # 分类表去前缀 | Remove prefixes from taxonomy table
tax = tax_names_repair(tax = tax, column_to_check = 7, column_to_add = 3)  # 分类表信息修复 | Repair taxonomy table names

# 数据预处理：数据抽平 | Data preprocessing: Rarefaction
otu_tax = data_rarefy(otu, method = "phyloseq", tax_table = tax)

# 拆分特征表和分类表 | Separate OTU table and taxonomy table
otu = otu_tax[["otu"]]
tax = otu_tax[["tax"]]

# 对齐代表性序列文件 | Align representative sequences file
otu_rep = merge(x = otu, y = rep, by = "#OTU ID", all.x = T, sort = F)

# 拆分特征表和代表性序列文件 | Separate OTU table and representative sequences
rep = otu_rep[, c(1, ncol(otu_rep))]
otu = otu_rep[, -c(ncol(otu_rep))]

# -----------------------------------------
# 数据分析与可视化 | Data Analysis & Visualization

# 物种堆叠图分析（门水平） | Stacked bar plot analysis (Phylum level)
data_sta_p = stackbar(otu, tax, metadata, tax_cla = "phylum", group1 = "group", group2 = "group2", row_n = 8)
# 可视化 | Visualization
stackbar_plot(data_sta_p, tax_cla = "phylum", title_legend = "Top 8 Phyla")

# 物种堆叠图分析（属水平） | Stacked bar plot analysis (Genus level)
data_sta_g = stackbar(otu, tax, metadata, tax_cla = "genus", group1 = "group", group2 = "group2", row_n = 20)
# 可视化 | Visualization
stackbar_plot(data_sta_g, tax_cla = "genus", title_legend = "Top 20 Genera")

# 弦图（门水平） | Chord diagram (Phylum level)
data_chord = chord(otu, metadata, tax, tax_cla = "phylum", group = "group2", row_n = 8)
# 可视化 | Visualization
chord_plot(data_chord)

# 韦恩图 | Venn diagram
data_venn = venn(otu, metadata, group = "group2")
# 可视化 | Visualization
venn_plot(data_venn)

# 集合图 | Upset plot
data_upset = Upset(otu, metadata, group = "group2")
# 可视化 | Visualization
Upset_plot(data_upset)

# 箱线图（Alpha 多样性分析） | Boxplot (Alpha diversity analysis)
data_alpha = alpha(otu, metadata, group = "group2", tree = tree)
# 可视化 | Visualization
alpha_plot(data_alpha)

# PCA | Principal Component Analysis (PCA)
data_pca = pca(otu, metadata, group = "group2")
# 可视化 | Visualization
pca_plot(data_pca)

# PCoA | Principal Coordinates Analysis (PCoA)
data_pcoa = pcoa(otu, metadata, group = "group2")
# 可视化 | Visualization
pcoa_plot(data_pcoa)

# NMDS | Non-metric Multidimensional Scaling (NMDS)
data_nmds = nmds(otu, metadata, group = "group2")
# 可视化 | Visualization
nmds_plot(data_nmds)

# RDA | Redundancy Analysis (RDA)
data_rda = RDA(otu, env, metadata, group = "group2")
# 可视化 | Visualization
RDA_plot(data_rda)

# CCA | Canonical Correspondence Analysis (CCA)
data_cca = CCA(otu, env, metadata, group = "group2")
# 可视化 | Visualization
CCA_plot(data_cca)

# 热图 | Heatmap
data_heatmap = heatmap(otu, tax, metadata, tax_cla = "genus", group1 = "group", group2 = "group2", row_n = 30)
# 可视化 | Visualization
heatmap_plot(data_heatmap, fontsize_col = 14, file_height = 10, file_width = 12)

# 共现性网络分析 | Co-occurrence network analysis
data_net = network(otu, tax, metadata, tax_cla = "genus")
data_net
```

Partial example figures:
> 部分示例图：

![image](https://github.com/user-attachments/assets/1a2e6ba6-6e1f-4846-97f7-81332c1d14b9)

<br>

## Additional Information 补充说明
### About the data_rarefy() function 
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
