# amplysis
![version](https://img.shields.io/badge/version-1.3.2-blue.svg)
![license](https://img.shields.io/badge/license-GPL--3.0-green.svg)
![platform](https://img.shields.io/badge/platform-R-75AADB.svg)
![status](https://img.shields.io/badge/status-active-brightgreen.svg)
![CRAN](https://img.shields.io/badge/CRAN-Submitting-orange.svg)

**An R package for rapid analysis of 16S amplicon sequencing data**

An R package for 16S rRNA gene amplicon sequencing data that integrates data preprocessing, analysis, and visualization methods. These methods include microbial composition analysis, α-diversity analysis, β-diversity analysis, differential analysis, correlation analysis, and network analysis, among others.

<br>

## Background

The rapid development of bioinformatics tools has enabled researchers to perform complex 16S rRNA gene amplicon sequencing analysis. However, for those without a bioinformatics background, existing tools and R packages can be complex and difficult to use. The amplysis project was created to address this issue, providing an accessible solution for researchers to easily conduct data analysis. 

<br>

## Install
### 1. Install the latest version from GitHub (recommended)

Install the R package `devtools`:

```
install.packages("devtools")
library(devtools)
```

Install the R package `amplysis` via `devtools`:

```
devtools::install_github("min-perilla/amplysis", upgrade = "never")
```

### 2. Install from a local `.tar.gz` source package

Go to the [Releases page](https://github.com/min-perilla/amplysis/releases) to download the latest `amplysis_*.tar.gz` release package and note its save location. Then use the following code in R to install the local source package (replace the path with your actual file location)


- Windows users can install by entering the following code in RStudio:

```
install.packages("C:/path/to/your/amplysis_X.X.X.tar.gz", repos = NULL, type = "source")
```

- macOS users can install by entering the following code in RStudio:
```
install.packages("/Users/path/to/your/amplysis_X.X.X.tar.gz", repos = NULL, type = "source")
```

### 3. Install from CRAN

> [!IMPORTANT]
> Submitting to CRAN, stay tuned...

```
install.packages("amplysis")
```

<br>

## Usage

Click [here](https://github.com/min-perilla/amplysis/releases/download/latest/example_data.zip) to download and extract the example_data.zip archive. After extraction, you will find the following files and folders:

* `test.R` (used to test the core functionality of the `amplysis` package)
* `test_readData.R` (used to test the data import function `read_data()`)
* `example data` folder (contains all types of example datasets)
* `other case studies data` folder (contains additional test datasets)

To test the `read_data()` function, open the `test_readData.R` file in `RStudio`, select all the code, and click `Run`. You will see how `read_data()` easily reads various types of data files.

To test the core functionality of `amplysis`, open the `test.R` file in `RStudio` select all the code, and click `Run`. You will see how the `amplysis` package quickly analyzes data and generates high-quality visualizations.

Alternatively, you can copy and paste the example code below into `RStudio` to try it out (Make sure the working directory is set properly):

```
# test.R
# Clear all variables
rm(list = ls())

# Load amplysis package
library(amplysis)

# Set working directory
set_wd()


# ------------------------------------------------------------------------------
# Load example data
otu = read_data("./example data/otu/otu.csv")            # OTU table
tax = read_data("./example data/tax/tax.csv")            # Taxonomy table
metadata = read_data("./example data/metadata/metadata.csv")  # Sample metadata
rep = read_data("./example data/rep_seqs/rep_seqs.csv")  # Representative sequences
env = read_data("./example data/env/env.csv")            # Environmental factors
tree = read_data("./example data/tree/tree_rooted.nwk")  # Phylogenetic tree


# ------------------------------------------------------------------------------
# Data preprocessing: Taxonomy table
tax = tax_separate(tax = tax, index = 2, delim = "; ")                     # Split taxonomy table into columns 
tax = tax_trim_prefix(tax = tax, index = c(2:8), length = 3)               # Split taxonomy table into columns
tax = tax_names_repair(tax = tax, column_to_check = 7, column_to_add = 3)  # Optimization of taxonomic annotation in the taxonomy table

# Data preprocessing: Rarefaction
otu_tax = data_rarefy(otu, method = "phyloseq", tax_table = tax)

# Separate OTU table and taxonomy table
otu = otu_tax[["otu"]]
tax = otu_tax[["tax"]]

# Align representative sequences file
otu_rep = merge(x = otu, y = rep, by = "#OTU ID", all.x = T, sort = F)

# Separate OTU table and representative sequences
rep = otu_rep[, c(1, ncol(otu_rep))]
otu = otu_rep[, -c(ncol(otu_rep))]


# ------------------------------------------------------------------------------
# Data Analysis & Visualization
# Stacked bar plot analysis (Phylum level)
data_sta_p = stackbar(otu, tax, metadata, tax_cla = "phylum", group1 = "group", group2 = "group2", row_n = 8)
stackbar_plot(data_sta_p, tax_cla = "phylum", title_legend = "Top 8 Phyla")  # Visualization

# Stacked bar plot analysis (Genus level)
data_sta_g = stackbar(otu, tax, metadata, tax_cla = "genus", group1 = "group", group2 = "group2", row_n = 20)
stackbar_plot(data_sta_g, tax_cla = "genus", title_legend = "Top 20 Genera")  # Visualization

# Chord diagram (Phylum level)
data_chord = chord(otu, metadata, tax, tax_cla = "phylum", group = "group2", row_n = 8)
chord_plot(data_chord)  # Visualization

# Venn diagram
data_venn = venn(otu, metadata, group = "group2")
venn_plot(data_venn)  # Visualization

# Upset plot
data_upset = Upset(otu, metadata, group = "group2")
Upset_plot(data_upset)  # Visualization

# Boxplot (Alpha diversity analysis)
data_alpha = alpha(otu, metadata, group = "group2", tree = tree)
alpha_plot(data_alpha)  # Visualization

# Principal Component Analysis (PCA)
data_pca = pca(otu, metadata, group = "group2")
pca_plot(data_pca)  # Visualization

# Principal Coordinates Analysis (PCoA)
data_pcoa = pcoa(otu, metadata, group = "group2")
pcoa_plot(data_pcoa)  # Visualization

# Non-metric Multidimensional Scaling (NMDS)
data_nmds = nmds(otu, metadata, group = "group2")
nmds_plot(data_nmds)  # Visualization

# Redundancy Analysis (RDA)
data_rda = RDA(otu, env, metadata, group = "group2")
RDA_plot(data_rda)  # Visualization

# Canonical Correspondence Analysis (CCA)
data_cca = CCA(otu, env, metadata, group = "group2")
CCA_plot(data_cca)  # Visualization

# Heatmap
data_heatmap = heatmap(otu, tax, metadata, tax_cla = "genus", group1 = "group", group2 = "group2", row_n = 30)
heatmap_plot(data_heatmap, fontsize_col = 14, file_height = 10, file_width = 12)  # Visualization

# Co-occurrence network analysis
# Visualization using Gephi is recommended
data_net = network(otu, tax, metadata, tax_cla = "genus")
data_net
```
<br>

Partial example figures:

![example](https://github.com/user-attachments/assets/9d722668-3e6f-4bd5-a1c3-894264c9cbd5)

<br>

## Additional Information
### 1. Preparing External Data for Analysis

If you already have downstream analysis files such as a feature table, taxonomy table, and representative sequences (rep-seqs), you can refer to the relevant content in this article (DOI link will be provided once available) and use the [example data files](https://github.com/min-perilla/amplysis/releases/download/latest/example_data.zip) as a reference. In most cases, only minimal or no modification is needed for direct use with the `amplysis` R package.


If you only have raw sequencing data (such as `fasta` or `fastq` files), we recommend visiting the [QIIME 2 official documentation](https://amplicon-docs.qiime2.org/en/latest/) to learn how to generate the downstream analysis data files through standard workflows.
<br>

### 2. Using `read_data()` to Read `.biom` Files

The `read_data()` function supports reading `JSON` format `.biom` files. If you already have a `.biom` file (such as one exported from `QIIME 2`, which may be in `HDF5` format), you can convert it to `JSON` format using the following methods.

Open the `Terminal` and enter the following command to check whether `biom` is installed:

```bash
biom --version
```

If not installed, use the following command to install the `biom-format` tool:

```bash
pip install biom-format
```

Convert the `HDF5` format `.biom` file to `JSON` format:

```bash
biom convert -i input_table.biom -o output_table_json.biom --to-json --table-type "OTU table"
```

If you attempt to read an `HDF5` format `.biom` file (such as one exported from `QIIME 2`), the `read_data()` function may produce the following warning:

```R
Warning messages:
1: In strsplit(conditionMessage(e), "\n") :
  unable to translate 'lexical error: invalid char in json text.
                                       <89>HDF                        (right here) ------^
' to a wide string
2: In strsplit(conditionMessage(e), "\n") : input string 1 is invalid
```
Although these warnings may appear, the file can still be read correctly. Rest assured, this issue will be resolved in a future update.
<br>

### 3. About the `data_rarefy()` function 

The R package provides a function `data_rarefy()` for data rarefaction. It integrates the `rrarefy()` function from the `vegan` package or the `rarefy_even_depth()` function from the `phyloseq` package for data rarefaction.  

Note that for the `rarefy_even_depth()` function in the `phyloseq` package, both replacement sampling and non-replacement sampling methods can be chosen. In replacement sampling, each time an element is selected from the sample pool, it is returned, keeping the sample pool unchanged, allowing some elements to be selected multiple times. This method can speed up computation and reduce memory usage, especially when handling large datasets. However, it may result in some OTU or ASV counts exceeding their original values, which could affect the accuracy of the analysis. 

In contrast, non-replacement sampling removes the selected sample after each draw, ensuring that each sample is selected only once, thereby ensuring OTU/ASV counts do not exceed their original values. The `rarefy_even_depth()` function from the `phyloseq` package enables replacement sampling by default (parameter `replace = TRUE`), which may have a slight impact on the resulting analysis. To ensure accuracy, the `data_rarefy()` function from the `amplysis` package uses non-replacement sampling for data rarefaction (parameter `replace = FALSE`).  

<br>

## Maintainers
[@min-perilla](https://github.com/min-perilla)

<br>

## Contributing
[@min-perilla](https://github.com/min-perilla)

<br>

## License
The project is licensed under **GPL (>= 3)**. For more details, please refer to the **[LICENSE](https://github.com/min-perilla/amplysis/blob/main/LICENSE)** file.

