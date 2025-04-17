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

# phylum
data_sta_p1 = stackbar(otu = otu, tax = tax, metadata = metadata, 
                      id_col = 1, tax_cla = "phylum", group1 = "group2", group2 = "group3", parallel_method = "none", row_n = 8)
stackbar_plot(data_sta_p1, tax_cla = "phylum", x_group = "group2", facet_group = "group3", 
              bar_width = 0.6, title_legend = "Top 8 Phyla", title_x = "Time",  
              filename = "phylum_1", file_height = 9, file_width = 16)

# genus
data_sta_g1 = stackbar(otu = otu, tax = tax, metadata = metadata, 
                      id_col = 1, tax_cla = "genus", group1 = "group2", group2 = "group3", parallel_method = "none", row_n = 20)
stackbar_plot(data_sta_g1, tax_cla = "genus", x_group = "group2", facet_group = "group3", 
              bar_width = 0.6, title_legend = "Top 20 Genera", title_x = "Time",  
              filename = "genus_1", file_height = 9, file_width = 16)

#------
# phylum
data_sta_p2 = stackbar(otu = otu, tax = tax, metadata = metadata, 
                      id_col = 1, tax_cla = "phylum", group1 = "group2", group2 = "group", parallel_method = "none", row_n = 8)
stackbar_plot(data_sta_p2, tax_cla = "phylum", x_group = "group2", facet_group = "group", bar_width = 0.6, title_legend = "Top 8 Phyla", 
              angle_x = 45, angle_hjust = 1, angle_vjust = 1, title_x = "Body Sites",  
              filename = "phylum_2", file_height = 9, file_width = 18)

# genus
data_sta_g2 = stackbar(otu = otu, tax = tax, metadata = metadata, 
                       id_col = 1, tax_cla = "genus", group1 = "group2", group2 = "group", parallel_method = "none", row_n = 20)
stackbar_plot(data_sta_g2, tax_cla = "genus", x_group = "group2", facet_group = "group", bar_width = 0.6, title_legend = "Top 20 Genera", 
              angle_x = 45, angle_hjust = 1, angle_vjust = 1, title_x = "Body Sites",  
              filename = "genus_2", file_height = 9, file_width = 18)

# PCoA
data_pcoa_1 = pcoa(otu, metadata, group = "group", parallel_method = "none")
pcoa_plot(data_pcoa_1, ellipse_type = 1)

