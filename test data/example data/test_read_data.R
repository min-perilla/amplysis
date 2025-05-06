# Clear all variables
rm(list = ls())

# Load R package
library(amplysis)

# Set working directory
set_wd()


## import data
# otu
otu_csv = read_data("./otu/otu.csv")
otu_tsv = read_data("./otu/otu.tsv")
otu_txt = read_data("./otu/otu.txt")
otu_xls = read_data("./otu/otu.xls")
otu_xlsx = read_data("./otu/otu.xlsx")
otu_biom_json = read_data("./otu/otu_json.biom")
otu_biom_hdf5 = read_data("./otu/otu_hdf5.biom")

# tax
tax_csv = read_data("./tax/tax.csv")
tax_tsv = read_data("./tax/tax.tsv")
tax_txt = read_data("./tax/tax.txt")
tax_xls = read_data("./tax/tax.xls")
tax_xlsx = read_data("./tax/tax.xlsx")

# rep_seqs
rep_csv = read_data("./rep_seqs/rep_seqs.csv")
rep_tsv = read_data("./rep_seqs/rep_seqs.tsv")
rep_txt = read_data("./rep_seqs/rep_seqs.txt")
rep_xls = read_data("./rep_seqs/rep_seqs.xls")
rep_xlsx = read_data("./rep_seqs/rep_seqs.xlsx")
rep_fasta = read_data("./rep_seqs/rep_seqs.fna")

# tree
tree = read_data("./tree/tree_rooted.nwk")

# env
env_csv = read_data("./env/env.csv")
env_tsv = read_data("./env/env.tsv")
env_txt = read_data("./env/env.txt")
env_xls = read_data("./env/env.xls")
env_xlsx = read_data("./env/env.xlsx")

# metadata
metadata_csv = read_data("./metadata/metadata.csv")
metadata_tsv = read_data("./metadata/metadata.tsv")
metadata_txt = read_data("./metadata/metadata.txt")
metadata_xls = read_data("./metadata/metadata.xls")
metadata_xlsx = read_data("./metadata/metadata.xlsx")


