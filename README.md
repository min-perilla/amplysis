# amplysis 1.3.0
**An R package for rapid analysis of 16S amplicon sequencing data**

> **16S 扩增子测序数据分析 R 包**

A series of analytical, statistical, and plotting methods in 16S rRNA gene amplicon sequencing data. These methods are used for data preprocessing, microbial composition analysis, alpha diversity analysis, beta diversity analysis, differential analysis, correlation and network analysis, among others.

> 一系列用于 16S rRNA 基因扩增子测序数据的分析、统计和绘图方法。这些方法用于数据预处理、微生物组成分析、α 多样性分析、β 多样性分析、差异分析、关联分析和网络分析等。

## Background  背景

The rapid development of bioinformatics tools has enabled researchers to perform complex 16S rRNA gene amplicon sequencing analysis. However, for those without a bioinformatics background, existing tools and R packages can be complex and difficult to use. The amplysis project was created to address this issue, providing an accessible solution for researchers to easily conduct data analysis. 

> 生物信息学工具的快速发展，使得研究者能够进行复杂的 16S rRNA 基因扩增子测序数据分析。然而，对于缺乏生物信息学背景的研究者来说，现有的工具和 R 包操作复杂，使用门槛较高。amplysis 集合 R 包项目应运而生，旨在为研究者提供一个易于使用的数据分析解决方案。

## Install  安装

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

```
install.packages("devtools")
devtools::install_github("min-perilla/amplysis")
```

## Usage  使用方法
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
chord_plot(data_chord)

# 韦恩图 | Venn diagram
data_venn = venn(otu, metadata, group = "group2")
venn_plot(data_venn)

# 集合图 | Upset plot
data_upset = Upset(otu, metadata = , group = "group")
Upset_plot(data_upset)

# 箱线图（Alpha 多样性分析） | Boxplot (Alpha diversity analysis)
data_alpha = alpha(otu, metadata, group = "group2", tree = tree)
alpha_plot(data_alpha)

# PCA | Principal Component Analysis (PCA)
data_pca = pca(otu, metadata, group = "group2")
pca_plot(data_pca)

# PCoA | Principal Coordinates Analysis (PCoA)
data_pcoa = pcoa(otu, metadata, group = "group2")
pcoa_plot(data_pcoa)

# NMDS | Non-metric Multidimensional Scaling (NMDS)
data_nmds = nmds(otu, metadata, group = "group2")
nmds_plot(data_nmds)

# RDA | Redundancy Analysis (RDA)
data_rda = RDA(otu, env, metadata, group = "group2")
RDA_plot(data_rda)

# CCA | Canonical Correspondence Analysis (CCA)
data_cca = CCA(otu, env, metadata, group = "group2")
CCA_plot(data_cca)

# 热图 | Heatmap
data_heatmap = heatmap(otu, tax, metadata, tax_cla = "genus", group1 = "group", group2 = "group2")
heatmap_plot(data_heatmap)

# 共现性网络分析 | Co-occurrence network analysis
data_net = network(otu, tax, metadata, tax_cla = "genus")
data_net
```


## Maintainers  项目主要负责人

[@min-perilla](https://github.com/min-perilla)

## Contributing  贡献
[@min-perilla](https://github.com/min-perilla)


## License  许可证

The project is licensed under **GPL (>= 3)**. For more details, please refer to the **[LICENSE.md](https://github.com/min-perilla/amplysis/blob/main/LICENSE.md)** file.

> 该项目采用 **GPL (>= 3)** 许可证，详情请参阅 **[LICENSE.md](https://github.com/min-perilla/amplysis/blob/main/LICENSE.md)** 文件。

