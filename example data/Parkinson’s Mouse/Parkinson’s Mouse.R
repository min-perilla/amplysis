rm(list = ls())

library(amplysis)

set_wd()

# import data
otu = read_data("table.tsv")
tax = read_data("taxonomy.tsv")
metadata = read_data("metadata.tsv")
tree = read_data("tree_rooted.nwk")

# tax processing
tax1 = tax_separate(tax, 2, "; ")
tax2 = tax_trim_prefix(tax1, c(2:8), 3)
tax3 = tax_names_repair(tax2, 7, 3)
tax3[is.na(tax3)] <- ""
tax3[tax3 == ""] <- "Unknown"
write.csv(tax3, "tax.csv", row.names = F)
tax = tax3

# Data Analysis and Visualization


# genus
data_sta_g1 = stackbar(otu = otu, tax = tax, metadata = metadata, 
                      id_col = 1, tax_cla = "genus", group1 = "group", group2 = "group2", parallel_method = "none", row_n = 20)
stackbar_plot(data_sta_g1, tax_cla = "genus", x_group = "group", facet_group = "group2", 
              bar_width = 0.6, title_legend = "Top 20 Genera", 
              filename = "genus_1", file_height = 9, file_width = 14)


# PCoA
data_pcoa_1 = pcoa(otu, metadata, group = "group", parallel_method = "none")
pcoa_plot(data_pcoa_1, ellipse_type = 1)

