#' @title R6 Class: amplysis
#'
#' @description
#' The R6 class `amplysis` encapsulates data preprocessing and data
#' analysis-related functions into a single class. This encapsulation simplifies
#' data management and preprocessing tasks during the analysis process. By using
#' the `amplysis` class, researchers can perform data input, modification,
#' preprocessing, and analysis operations through a unified interface, thereby
#' enhancing work efficiency.
#'
#' @return R6 Class: otu, tax, rep, tree, env, metadata, data, chart
#'
#' @export
#'
#' @importFrom R6 R6Class
#'
#'
amplysis = R6::R6Class(
  # Define an R6 class: Base class.
  "amplysis",


  ##
  # inherit
  inherit = otu_tax_rep_tree_env_metadata_data_chart,


  ##
  public = list(
    #' @description
    #' initialization function.
    #'
    #' @param otu otu table
    #' @param tax tax table
    #' @param rep_seq representative sequence
    #' @param tree phylogenetic tree file
    #' @param env environmental factor
    #' @param metadata metadata table
    #'
    #' @return R6 class
    #'
    #' @examples \dontrun{amplysis$new(otu, tax, rep_seq, tree, env, metadata)}
    #'
    initialize = function(otu = NULL, tax = NULL, rep_seq = NULL, tree = NULL, env = NULL, metadata = NULL) {
      private$.otu = otu
      private$.tax = tax
      private$.rep_seq = rep_seq
      private$.tree = tree
      private$.env = env
      private$.metadata = metadata
      private$.data = NULL
      private$.chart = NULL
    },


    #' @description
    #' verride the print function.
    #'
    #' @param ...
    #'
    #' @return gg, ggplot2, chart
    #'
    #' @examples \dontrun{
    #' a <- amplysis$new(otu, metadata)
    #' print(a)}
    #'
    print = function(...) {
      print(private$.chart)
    },


    ## Data preprocessing
    ## tax_separate
    #' @description
    #' Call the `separate_wider_delim()` function from the tidyr package to
    #' separate the data in the taxonomy column of the data table.
    #'
    #' @param index (integer) The column number of the taxonomy column in the data table.
    #' @param delim (character) The delimiter. Default is ";".
    #' @param names (character) Names for each column after data separation.
    #' Ensure the number of names matches the number of columns after separation,
    #' such as 'names = c("domain", "phylum", "class", "order", "family", "genus", "species")'.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_separate(index = 2)}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_separate(index = 2, delim = ";", names = c("domain", "phylum",
    #' "class", "order", "family", "genus", "species"))}
    #'
    tax_separate = function(index, delim, names) {
      private$.tax <- amplysis::tax_separate(
        tax = private$.tax,
        index = index,
        delim = delim,
        names = names
      )
    },


    ## tax_trim_prefix
    #' @description
    #' The function can remove a fixed-length prefix from taxonomic information, such as p__Firmicutes -> Firmicutes.
    #'
    #' @param index (integer or character) Which columns contain the classification information to be repaired,
    #' for example: index = 3, index = c(2:8), index = c("domain", "phylum", "class", "order", "family", "genus", "species").
    #' @param position (integer) The length of the prefix.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_trim_prefix(index = c(2:8), position = 3)}
    #'
    tax_trim_prefix = function(index, position) {
      private$.tax = amplysis::tax_trim_prefix(
        tax = private$.tax,
        index = index,
        position = position
      )
    },


    ## tax_names_repair
    #' @description
    #' In taxonomic information, terms such as \"unclassified\", \"unknown\", \"uncultured\",
    #' \"norank\", \"unidentified\", \"Unknown_Species\", \"Unknown_Genus\", \"Unknown_Family\",
    #' \"Unknown_Order\", \"Unknown_Class\", and \"Unknown_Phylum\" may exist.
    #' When plotting a species composition chart, simply merging these categories may result in the loss of some taxonomic information.
    #'
    #' @description
    #' For instance, one entry may be categorized as Firmicutes at the phylum level and as
    #' \"uncultured\" at the genus level, while another entry may be categorized as Armatimonadetes at
    #' the phylum level and also as \"uncultured\" at the genus level. Merging them at
    #' the genus level would result in both entries being classified as \"uncultured\",
    #' leading to the loss of some taxonomic information. The correct approach is to refine
    #' and repair the taxonomic information. In this case, the entries at the genus level should be
    #' specified as \"uncultured_Firmicutes\" and \"uncultured_Armatimonadetes\", respectively,
    #' to preserve the detailed taxonomic information.
    #'
    #' @param column_to_check (integer or character) The columns to be repaired, such as column_to_check = 8, column_to_check = c(4:8),
    #' If column names are known, you can also input column_to_check = "column_name"
    #' @param column_to_add (integer or character) If column_to_check is unclassified, append the column number as a suffix.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_names_repair(column_to_check = 7, column_to_add = 3)}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_names_repair(column_to_check = c(4:8), column_to_add = 3)}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_names_repair(column_to_check = "genus", column_to_add = "phylum")}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_names_repair(column_to_check = c("genus", "species"), column_to_add = "phylum")}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$tax_names_repair(column_to_check = c("genus", "species"), column_to_add = 3)}
    #'
    tax_names_repair = function(column_to_check, column_to_add) {
      private$.tax = amplysis::tax_names_repair(
        tax = private$.tax,
        column_to_check = column_to_check,
        column_to_add = column_to_add)
    },


    ## align_otu_tax
    #' @description
    #' When using the rarefaction method of the phyloseq package with the parameter
    #' trimOTUs = TRUE, some OTUs in the OTU table will be removed. In this case,
    #' the number of rows in the OTU table and the TAX table will not be consistent.
    #' The align_otu_tax() function can be used to align these two tables.
    #'
    #' @param by.x (integer or character) Used to specify the column name in the tax data frame,
    #' defaulting to 1, indicating the column name of the first column. Alternatively,
    #' you can directly enter the column name, such as by.x = "OTU ID".
    #' @param by.y (integer or character) Used to specify the column name in the otu data frame,
    #' defaulting to 1, indicating the column name of the first column. Alternatively,
    #' you can directly enter the column name, such as by.x = "OTU ID".
    #' @param all.x (logical) Should all rows of the tax data frame be retained? Default is FALSE.
    #' @param all.y (logical) Should all rows of the otu data frame be retained? Default is FALSE.
    #' @param sort (logical) Should reordering be done? Default is no.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$align_otu_tax(by.x = 1, by.y = 1, all.x = F, all.y = T, sort = F)}
    #'
    align_otu_tax = function(by.x = 1, by.y = 1, all.x = F, all.y = F, sort = F) {
      otu_tax = amplysis::align_otu_tax(
        tax_table = private$.tax,
        otu_table = private$.otu,
        by.x = by.x,
        by.y = by.y,
        all.x = all.x,
        all.y = all.y,
        sort = sort
      )

      #
      private$.otu = otu_tax[["otu"]]
      private$.tax = otu_tax[["tax"]]
    },


    ## data_rarefy
    #' @description
    #' Use the "rrarefy" function from the "vegan" package or the "rarefy_even_depth" function from the "phyloseq" package for data rarefaction.
    #'
    #' @param id_col (integer) The OTU_ID column is in which column, defaulting to 0 means there is no OTU_ID column, and the data is already numeric.
    #' @param method (character) The data rarefaction can be performed using the "vegan" package or the "phyloseq" package.
    #' @param seed (integer) (Only applicable when using the "phyloseq" package for data rarefaction) Random seed.
    #' @param replace (logical) Whether to sample with replacement (TRUE) or without replacement (FALSE).
    #' (1) Sampling with replacement is faster and more memory-efficient;
    #' (2) Sampling with replacement means that the read counts for some OTUs could
    #' be greater than the original count values.
    #' @param trimOTUs (logical) (Only applicable when using the "phyloseq" package for data rarefaction) Whether to remove OTUs with zero counts for each sample.
    #' @param tax_table (tax table)
    #' @param write_file (logical) Whether to save the rarefied file.
    #' @param file_name (character) Customize the saved file name.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$data_rarefy()}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$data_rarefy(id_col = 0, method = "vegan")}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$data_rarefy(id_col = 0, method = "phyloseq", seed = 1231)}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$data_rarefy(id_col = 1, method = "phyloseq", trimOTUs = T,
    #' tax_table = tax, write_file = T)}
    #'
    data_rarefy = function(id_col = 0, method = "vegan", seed = 123, replace = F,
                           trimOTUs = T, tax_table = NULL, write_file = TRUE,
                           file_name = "table"){
      result = amplysis::data_rarefy(
        file = private$.otu,
        id_col = id_col,
        method = method,
        seed = seed,
        replace = replace,
        trimOTUs = trimOTUs,
        tax_table = tax_table,
        write_file = write_file,
        file_name = file_name
      )

      ##
      if(is.null(tax_table)){
        private$.otu = result
      } else {
        private$.otu = result[[1]]
        private$.tax = result[[2]]
      }
    },



    ## data analysis
    ## alpha
    #' @description
    #' Alpha diversity analysis includes the calculation of Shannon, Simpson, Chao1,
    #' Ace, Pielou, Goods_coverage, and PD indices. It supports significance
    #' analysis using one-way ANOVA and multiple comparison methods.
    #' The significance marking method uses the letter marking method.
    #'
    #' @param id_col (integer) The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param group (character) Group 1, please enter the column name of
    #' the grouping information in the metadata table.
    #' @param parallel_method (character) Sample processing methods for the same group:
    #' average, sum, median, none.
    #' @param method Multiple comparison methods, please enter a number from 1 to 7
    #' to choose:
    #'
    #' "1" = "Tukey-HSD test",
    #'
    #' "2" = "Fisher-LSD test",
    #'
    #' "3" = "Student-Newman-Keuls test",
    #'
    #' "4" = "Duncan test",
    #'
    #' "5" = "Scheffe test",
    #'
    #' "6" = "Waller-Duncan test",
    #'
    #' "7" = "REGW test".
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$alpha(id_col = 1, group = "group", parallel_method = "mean",
    #' tree = tree, method = 1)
    #' }
    #'
    alpha = function(id_col = 1, group = "group",
                     parallel_method = "mean", method = 1)
    {
      result = amplysis::alpha(
        otu = private$.otu,
        metadata = private$.metadata,
        id_col = id_col,
        group = group,
        parallel_method = parallel_method,
        tree = private$.tree,
        method = method
      )


      ##
      private$.data = result

    },


    ## alpha_plot
    #' @description
    #' After processing with the `alpha()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #' @param custom_order (character) Custom legend order.
    #'
    #' @param size_point (numeric) The size of points.
    #' @param size_differ (numeric) The size of significance marker letters.
    #' @param errorbar_width (numeric) The width of the horizontal lines on the error bars.
    #' @param errorbar_linewidth (numeric) The width of the vertical lines on the error bars.
    #'
    #' @param title (character) The title of the chart.
    #' @param title_x (character) The title of the X-axis.
    #' @param title_y (character) The title of the Y-axis.
    #'
    #' @param size_title (numeric) Font size of the main title.
    #' @param size_title_x (numeric) Font size of the horizontal axis title.
    #' @param size_title_y (numeric) Font size of the vertical axis title.
    #'
    #' @param size_x (numeric) Font size of the horizontal axis tick labels.
    #' @param size_y (numeric) Font size of the vertical axis tick labels.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$alpha(id_col = 1, group = "group", parallel_method = "mean", tree = tree, method = 1)
    #' a$alpha_plot(data,
    #' color_scheme = c('#aaf200','#0082ff',"#d23aa4","#c777ff", "#79ff79"),
    #' custom_order = c("A", "B", "R", "D", "S"))
    #' }
    #'
    alpha_plot = function(color_scheme = NULL, custom_order = NULL,
                          size_point = 5, size_differ = 14, errorbar_width = 0.15,
                          errorbar_linewidth = 0.8,
                          title = NULL, title_x = NULL, title_y = NULL,
                          size_title = 40, size_title_x = 28, size_title_y = 28,
                          size_x = 28, size_y = 28,
                          filename = "alpha", file_width = 12, file_height = 9)
    {
      private$.chart = amplysis::alpha_plot(
        data = private$.data, color_scheme = color_scheme, custom_order = custom_order,
        size_point = size_point, size_differ = size_differ, errorbar_width = errorbar_width,
        errorbar_linewidth = errorbar_linewidth,
        title = title, title_x = title_x, title_y = title_y,
        size_title = size_title, size_title_x = size_title_x, size_title_y = size_title_y,
        size_x = size_x, size_y = size_y,
        filename = filename, file_width = file_width, file_height = file_height
      )
      print(private$.chart)
    },
    ############################################################################


    ## RDA
    #' @description
    #' Redundancy Analysis (RDA) is a statistical method used to explore the
    #' relationship between response variables and explanatory variables in
    #' multivariate datasets. This function references the `rda()` function in the
    #' `vegan` package for analysis. Prior to conducting RDA analysis, it is
    #' recommended to perform Detrended Correspondence Analysis (DCA) detection
    #' using the `decorana` function from the `vegan` package. Selection is based
    #' on the Axis Lengths value of DCA1 as follows: (1) If >4.0, choose CCA;
    #' (2) If between 3.0-4.0, both RDA and CCA are acceptable;
    #' (3) If <3.0, choose RDA.
    #'
    #' @param id_col (integer) The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param group (character) Group 1, please enter the column name of
    #' the grouping information in the metadata table.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$RDA(1, "group")
    #' }
    #'
    RDA = function(id_col = 1, group = "group") {
      private$.data = amplysis::RDA(
        otu = private$.otu,
        env = private$.env,
        metadata = private$.metadata,
        id_col = id_col,
        group = group
      )
    },


    ## RDA_plot
    #' @description
    #' After processing with the `RDA()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #' @param group (character) Grouping information.
    #' @param custom_order Custom legend order.
    #' @param seed Seed
    #'
    #' @param size_point (numeric) The size of points.
    #' @param size_point_legend (numeric) The size of legend points.
    #' @param spacing_legend_point (numeric) The internal padding of the legend.
    #' @param spacing_legend_title (numeric) The spacing between the legend title and the body
    #' @param legend_ncol (integer) Number of columns in the legend.
    #'
    #' @param title (character) Main title.
    #' @param title_legend (character) Legend title.
    #'
    #' @param size_title (numeric) Font size of the main title.
    #' @param size_title_x (numeric) Font size of the horizontal axis title.
    #' @param size_title_y (numeric) Font size of the vertical axis title.
    #' @param size_title_legned (numeric) Font size of legend title.
    #'
    #' @param size_x (numeric) Font size of the horizontal axis tick labels.
    #' @param size_y (numeric) Font size of the vertical axis tick labels.
    #' @param size_legned (numeric) Font size of the legend.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$RDA(1, "group")
    #' a$RDA_plot(group = "group", custom_order = NULL,
    #'            color_scheme = c("#00b0f6", "#FFC24B", "#f8766d"), seed = 123,
    #'
    #'            size_point = 8, size_point_legend = 8, spacing_legend_point = 1.2,
    #'            spacing_legend_title = 0.5, legend_ncol = 1,
    #'
    #'            title = "RDA", title_legend = "Group",
    #'
    #'            size_title = 28, size_title_x = 20, size_title_y = 20,
    #'            size_title_legned = 24,
    #'
    #'            size_x = 16, size_y = 16, size_legned = 16,
    #'
    #'            filename = "RDA", file_width = 16, file_height = 9)
    #' }
    #'
    RDA_plot = function(color_scheme = NULL, group = "group", custom_order = NULL, seed = 123,

                        size_point = 4.5, size_point_legend = 8, spacing_legend_point = 1.2,
                        spacing_legend_title = 0.5, legend_ncol = 1,

                        title = "RDA", title_legend = "Group",

                        size_title = 28, size_title_x = 20, size_title_y = 20,
                        size_title_legned = 24,

                        size_x = 18, size_y = 18, size_legned = 16,

                        filename = "RDA", file_width = 16, file_height = 9)
    {
      private$.chart = amplysis::RDA_plot(
        data = private$.data,
        color_scheme = color_scheme,
        group = group,
        custom_order = custom_order,
        seed = seed,

        size_point = size_point,
        size_point_legend = size_point_legend,
        spacing_legend_point = spacing_legend_point,

        spacing_legend_title = spacing_legend_title,
        legend_ncol = legend_ncol,


        title = title,
        title_legend = title_legend,


        size_title = size_title,
        size_title_x = size_title_x,
        size_title_y = size_title_y,

        size_title_legned = size_title_legned,


        size_x = size_x,
        size_y = size_y,
        size_legned = size_legned,


        filename = filename,
        file_width = file_width,
        file_height = file_height)
      print(private$.chart)
    },
    ############################################################################


    ## CCA
    #' @description
    #' Canonical Correspondence Analysis (CCA) is a multivariate statistical method
    #' used to study the relationship between two sets of variables. It can reveal
    #' the inherent connections between the two sets of variables. This function
    #' references the `cca()` function in the `vegan` package for analysis. Prior to
    #' conducting RDA analysis, it is recommended to perform Detrended
    #' Correspondence Analysis (DCA) detection using the `decorana` function from
    #' the `vegan` package. Selection is based on the Axis Lengths value of DCA1 as
    #' follows: (1) If >4.0, choose CCA; (2) If between 3.0-4.0, both RDA and CCA
    #' are acceptable; (3) If <3.0, choose RDA.
    #'
    #'
    #' @param id_col (integer)The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param group (character) Group 1, please enter the column name of
    #' the grouping information in the metadata table.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$CCA(metadata, 1, "group")
    #' }
    #'
    CCA = function(id_col = 1, group = "group") {
      private$.data = amplysis::CCA(
        otu = private$.otu,
        env = private$.env,
        metadata = private$.metadata,
        id_col = id_col,
        group = group
      )
    },
    ############################################################################


    ## CCA_plot
    #' @description
    #' After processing with the `CCA()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #' @param group (character) Grouping information.
    #' @param custom_order Custom legend order.
    #' @param seed Seed
    #'
    #' @param size_point (numeric) The size of points.
    #' @param size_point_legend (numeric) The size of legend points.
    #' @param spacing_legend_point (numeric) The internal padding of the legend.
    #' @param spacing_legend_title (numeric) The spacing between the legend title and the body
    #' @param legend_ncol (integer) Number of columns in the legend.
    #'
    #' @param title (character) Main title.
    #' @param title_legend (character) Legend title.
    #'
    #' @param size_title (numeric) Font size of the main title.
    #' @param size_title_x (numeric) Font size of the horizontal axis title.
    #' @param size_title_y (numeric) Font size of the vertical axis title.
    #' @param size_title_legned (numeric) Font size of legend title.
    #'
    #' @param size_x (numeric) Font size of the horizontal axis tick labels.
    #' @param size_y (numeric) Font size of the vertical axis tick labels.
    #' @param size_legned (numeric) Font size of the legend.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$CCA(id_col = 1, group = "group")
    #' a$CCA_plot(group = "group", custom_order = NULL,
    #'            color_scheme = c("#00b0f6", "#FFC24B", "#f8766d"), seed = 123,
    #'
    #'            size_point = 8, size_point_legend = 8, spacing_legend_point = 1.2,
    #'            spacing_legend_title = 0.5, legend_ncol = 1,
    #'
    #'            title = "CCA", title_legend = "Group",
    #'
    #'            size_title = 28, size_title_x = 20, size_title_y = 20,
    #'            size_title_legned = 24,
    #'
    #'            size_x = 16, size_y = 16, size_legned = 16,
    #'
    #'            filename = "CCA", file_width = 16, file_height = 9)
    #' }
    #'
    CCA_plot = function(color_scheme = NULL, group = "group", custom_order = NULL, seed = 123,
                        size_point = 4.5, size_point_legend = 8, spacing_legend_point = 1.2,
                        spacing_legend_title = 0.5, legend_ncol = 1,
                        title = "RDA", title_legend = "Group",
                        size_title = 28, size_title_x = 20, size_title_y = 20,
                        size_title_legned = 24,
                        size_x = 18, size_y = 18, size_legned = 16,
                        filename = "RDA", file_width = 16, file_height = 9)
    {
      private$.chart = amplysis::CCA_plot(
        data = private$.data,
        color_scheme = color_scheme, group = group, custom_order = custom_order, seed = seed,

        size_point = size_point, size_point_legend = size_point_legend, spacing_legend_point = spacing_legend_point,
        spacing_legend_title = spacing_legend_title, legend_ncol = legend_ncol,

        title = title, title_legend = title_legend,

        size_title = size_title, size_title_x = size_title_x, size_title_y = size_title_y,
        size_title_legned = size_title_legned,

        size_x = size_x, size_y = size_y, size_legned = size_legned,
        filename = filename, file_width = file_width, file_height = file_height
        )
      print(private$.chart)
    },
    ############################################################################


    ## chord
    #' @description
    #' Process the feature table into the necessary data format for drawing a chord
    #' diagram.
    #'
    #' @param id_col (integer) The OTU_ID column is in which column,
    #' defaulting to 0 means there is no OTU_ID column, and the data is already
    #' numeric.
    #' @param tax_cla (character) Taxonomic level. Only column names from the tax
    #' table can be entered, for example tax_cla = 'genus'.
    #'
    #' @param group (Required, character) Grouping information. please enter the column name of
    #' the grouping information in the metadata table.
    #' @param parallel_method  (character) Sample processing methods for the same group:
    #' average, sum, median, none.
    #' @param row_n (integer) Preserve the top N taxa (including the Nth) based on
    #' abundance, while merging taxa with lower abundance into "others".
    #' @param digits (integer) When the parallel sample method is set to "mean," the
    #' number of decimal places will default to 0.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$chord(id_col = 1, tax_cla = "genus", group = "group",
    #' parallel_method = "mean", row_n = 8, digits = 0)}
    #'
    chord = function(id_col = 1, tax_cla = "genus",
                     group = "group", parallel_method = "mean", row_n = 8, digits = 0){
      private$.data = amplysis::chord(
        otu = private$.otu,
        metadata = private$.metadata,
        tax = private$.tax,
        id_col = id_col,
        tax_cla = tax_cla,
        group = group,
        parallel_method = parallel_method,
        row_n = row_n,
        digits = digits
      )
    },


    ## chord_plot
    #' @description
    #' After processing with the `chord()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #'
    #' @param size_axis (numeric) Font size on the sector axis.
    #' @param size_label (numeric) Font size of sector labels.
    #' @param label_height (numeric) Track height of label text.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$chord(id_col = 1, tax_cla = "genus", group = "group",
    #' parallel_method = "mean", row_n = 8, digits = 0)
    #' a$chord_plot(color_scheme = c("#27e6ff", "#42ff0e", "#33BEB7", "#F66320",
    #' "#FBA127", "#A463D7", "#DB3937", "#ffaec8", "#828282"))
    #' }
    #'
    chord_plot = function(color_scheme = NULL, size_axis = 0.8,
                          size_label = 1, label_height = 0.35, filename = "chord",
                          file_width = 16, file_height = 10)
    {
      private$.chart = amplysis::chord_plot(
        data = private$.data,
        color_scheme = color_scheme,
        size_axis = size_axis,
        size_label = size_label,
        label_height = label_height,
        filename = filename,
        file_width = file_width,
        file_height = file_height
      )
      print(private$.chart)
    },
    ############################################################################


    ## heat_map
    ## heatmap
    #' @description
    #' To perform heatmap analysis, you will need three files:
    #' 1. OTU table (containing abundance information);
    #' 2. Taxonomy table (containing species annotation information);
    #' 3. Metadata table (containing grouping information).
    #' This function undergoes a series of processing steps,
    #' such as clustering based on specified taxonomic levels using column numbers
    #' or names from the tax table, to generate a dataframe directly usable for
    #' pheatmap plotting. For quick plotting, you can use the built-in
    #' `heat_map_plot()` function to generate high-quality plots efficiently.
    #'
    #' @param id_col (integer) The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param tax_cla (character) Taxonomic level. Only column names from the tax
    #' table can be entered, for example tax_cla = 'genus'.
    #' @param group (character) Group 1, please enter the column name or column
    #' number of the grouping information in the metadata table.
    #' @param group2 (character) Group 2 for facetting plots, please enter the
    #' column name or column number of the grouping information in the metadata
    #' table.
    #' @param parallel_method (character) Parallel sample processing method,
    #' defaulting to mean. Options: mean (average), sum (summation),
    #' median (median).
    #' @param row_n (integer) Preserve the top N taxa (including the Nth) based on
    #' abundance.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$heat_map(id_col = 1, group = "group", group2 = NULL, tax_cla = "genus",
    #' parallel_method = "mean", row_n = 50)}
    #'
    heat_map = function(id_col = 1, tax_cla = "genus", group = "group",
                        group2 = NULL, parallel_method = "mean", row_n = 35)
    {
      private$.data = amplysis::heat_map(
        otu = private$.otu,
        tax = private$.tax,
        metadata = private$.metadata,
        id_col = id_col,
        tax_cla = tax_cla,
        group = group,
        group2 = group2,
        parallel_method = parallel_method,
        row_n = row_n
      )
    },


    ## heat_map_plot
    #' @description
    #' After performing heatmap analysis using the `heat_map()` function,
    #' you will obtain a dataframe. This function allows you to visualize the
    #' dataframe, thereby generating a high-quality heatmap.
    #' If you want to customize the x-axis order of the heatmap, please input a
    #' vector into the parameter cluster_cols. For example,
    #' cluster_cols = c(4, 5, 6, 1, 2, 3), indicates sorting according to column
    #' numbers 4 5 6 1 2 3.
    #'
    #' @param scale (character) scale is used to set normalization. "row" represents
    #' row-wise normalization, "column" represents column-wise normalization, and
    #' "none" represents no normalization.
    #' @param cellwidth (numeric) Translation: Represents the width of a single
    #' cell, default is "NA".
    #' @param cellheight (numeric) Translation: Represents the height of a single
    #' cell, default is "NA".
    #' @param color (character) The heatmap cell colors are generated automatically
    #' with a gradient.
    #' For example: c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24")
    #' @param cluster_cols (logical or integer vector) Whether to enable column
    #' clustering (when clustering is enabled,
    #' custom sorting is not possible). If you want to customize the x-axis order of
    #' the heatmap,
    #' please input a vector into the parameter cluster_cols.
    #' For example, cluster_cols = c(4, 5, 6, 1, 2, 3), indicates sorting according
    #' to column numbers 4 5 6 1 2 3.
    #' @param clustering_method (character) "Represents clustering methods,
    #' including: 'ward.D', ward.D2', 'single', 'complete', 'average', 'mcquitty',
    #' 'median', 'centroid'"
    #' @param treeheight_row (numeric) Row clustering tree height
    #' @param treeheight_col (numeric) Col clustering tree height
    #'
    #' @param filename (character) File name for saving
    #' @param file_width (numeric) Image width
    #' @param file_height (numeric) Image height
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$heat_map(id_col = 1, group = "group", group2 = NULL, tax_cla = "genus",
    #' parallel_method = "mean", row_n = 50)
    #' a$heat_map_plot(scale = "row", cellwidth = NA, cellheight = NA,
    #' color =  c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24"),
    #' cluster_cols = F, clustering_method = "ward.D", treeheight_row = 100,
    #' treeheight_col = 50, filename = "heatmap.png")}
    #'
    heat_map_plot = function(scale = "row", cellwidth = NA, cellheight = NA,
                             color =  c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24"),
                             cluster_cols = F, clustering_method = "ward.D", treeheight_row = 50,
                             treeheight_col = 50, filename = "heatmap.png", file_width = 12,
                             file_height = 12)
    {
      private$.chart = amplysis::heat_map_plot(
        mat = private$.data,
        scale = scale,
        cellwidth = cellwidth,
        cellheight = cellheight,
        color =  color,
        cluster_cols = cluster_cols,
        clustering_method = clustering_method,
        treeheight_row = treeheight_row,
        treeheight_col = treeheight_col,
        filename = filename,
        file_width = file_width,
        file_height = file_height
      )
      print(private$.chart)
    },
    ############################################################################


    ## pca
    #' @description
    #' In bioinformatics analysis, Principal Component Analysis (PCA) is a commonly
    #' used statistical technique for exploring and visualizing patterns and
    #' structures within high-dimensional datasets. The "pca()" function utilizes
    #' the built-in R function "stats::prcomp()" to analyze and calculate the
    #' explanatory power of each principal component.
    #'
    #' @param id_col (integer) The OTU_ID column is in which column,
    #' defaulting to 0 means there is no OTU_ID column, and the data is already
    #' numeric.
    #' @param group (Required, character) Grouping information. please enter the column name of
    #' the grouping information in the metadata table.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$pca(1, metadata)
    #' }
    #'
    pca = function(id_col = 1, group = "group")
    {
      private$.data = amplysis::pca(
        otu = private$.otu,
        metadata = private$.metadata,
        id_col = 1,
        group = "group"
      )
    },


    ## pca_plot
    #' @description
    #' After processing with the `pca()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #' @param custom_order (character) Custom legend order.
    #' @param seed Seed
    #'
    #' @param size_point (numeric) The size of points.
    #' @param size_point_legend (numeric) The size of legend points.
    #' @param spacing_legend_point (numeric) The internal padding of the legend.
    #' @param spacing_legend_title (numeric) The spacing between the legend title and the body
    #' @param legend_ncol (integer) Number of columns in the legend.
    #' @param label_is (logical) Whether to display data labels.
    #' @param size_label (numeric) Font size of the data label.
    #'
    #' @param title (character) Main title.
    #' @param title_sub (character) Subtitle.
    #' @param title_legend (character) Legend title.
    #'
    #' @param size_title (numeric) Font size of the main title.
    #' @param size_title_sub (numeric) Font size of the subtitle.
    #' @param size_title_x (numeric) Font size of the horizontal axis title.
    #' @param size_title_y (numeric) Font size of the vertical axis title.
    #' @param size_title_legned (numeric) Font size of legend title.
    #'
    #' @param size_x (numeric) Font size of the horizontal axis tick labels.
    #' @param size_y (numeric) Font size of the vertical axis tick labels.
    #' @param size_legned (numeric) Font size of the legend.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$pca(1, metadata)
    #' a$pca_plot(color_scheme = c("#00b0f6", "#FFC24B", "#f8766d", "#ae876d", "#AFC24B"),
    #'            seed = 123, custom_order = c("A", "B", "R", "D", "S"),
    #'
    #'            size_point = 8, size_point_legend = 8, spacing_legend_point = 1.2,
    #'            spacing_legend_title = 0.5, legend_ncol = 1,
    #'
    #'            title = "PCA", title_legend = "Group",
    #'
    #'            size_title = 28, size_title_x = 20, size_title_y = 20,
    #'            size_title_legned = 24,
    #'
    #'            size_x = 18, size_y = 18, size_legned = 16, filename = "PCA",
    #'            file_width = 12, file_height = 9)
    #' }
    #'
    #'
    pca_plot = function(color_scheme = NULL, custom_order = NULL, seed = 123,
                        size_point = 4.5, size_point_legend = 8, spacing_legend_point = 1.2,
                        spacing_legend_title = 0.5, legend_ncol = 1, label_is = T,
                        size_label = 5, title = "PCA", title_sub = NULL, title_legend = "Group",
                        size_title = 28, size_title_sub = 16, size_title_x = 20,
                        size_title_y = 20, size_title_legned = 24,
                        size_x = 18, size_y = 18, size_legned = 16,
                        filename = "PCA", file_width = 12, file_height = 9)
    {
      private$.chart = amplysis::pca_plot(
        data = private$.data,
        color_scheme = color_scheme,
        custom_order = custom_order,
        seed = seed,

        size_point = size_point,
        size_point_legend = size_point_legend,
        spacing_legend_point = spacing_legend_point,
        spacing_legend_title = spacing_legend_title,
        legend_ncol = legend_ncol,
        label_is = label_is,
        size_label = size_label,

        title = title,
        title_sub = title_sub,
        title_legend = title_legend,

        size_title = size_title,
        size_title_sub = size_title_sub,
        size_title_x = size_title_x,
        size_title_y = size_title_y,
        size_title_legned = size_title_legned,

        size_x = size_x,
        size_y = size_y,
        size_legned = size_legned,

        filename = filename,
        file_width = file_width,
        file_height = file_height
      )
      print(private$.chart)
    },
    ############################################################################


    ## pcoa
    #' @description
    #' Principal Coordinates Analysis(PCoA) is a visualization method for studying
    #' data similarities or dissimilarities, enabling the observation of differences
    #' between individuals or groups. The "pcoa()" function calculates the
    #' Bray-Curtis distance using "vegan::vegdist()", then performs PCoA using the
    #' built-in R function "stats::cmdscale()", etaining the eigenvalues.
    #'
    #' @param id_col (integer) The OTU_ID column is in which column,
    #' defaulting to 0 means there is no OTU_ID column, and the data is already
    #' numeric.
    #' @param group (Required, character) Grouping information. please enter the column name of
    #' the grouping information in the metadata table.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$pcoa(1, metadata)
    #' }
    #'
    pcoa = function(id_col = 1, group = "group")
    {
      private$.data = amplysis::pcoa(
        otu = private$.otu,
        metadata = private$.metadata,
        id_col = id_col,
        group = group
      )
    },


    ## pcoa_plot
    #' @description
    #' After processing with the `pcoa()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #' @param custom_order Custom legend order.
    #' @param seed Seed
    #'
    #' @param size_point (numeric) The size of points.
    #' @param size_point_legend (numeric) The size of legend points.
    #' @param spacing_legend_point (numeric) The internal padding of the legend.
    #' @param spacing_legend_title (numeric) The spacing between the legend title and the body
    #' @param legend_ncol (integer) Number of columns in the legend.
    #' @param label_is (logical) Whether to display data labels.
    #' @param size_label (numeric) Font size of the data label.
    #'
    #' @param title (character) Main title.
    #' @param title_sub (character) Subtitle.
    #' @param title_legend (character) Legend title.
    #'
    #' @param size_title (numeric) Font size of the main title.
    #' @param size_title_sub (numeric) Font size of the subtitle.
    #' @param size_title_x (numeric) Font size of the horizontal axis title.
    #' @param size_title_y (numeric) Font size of the vertical axis title.
    #' @param size_title_legned (numeric) Font size of legend title.
    #'
    #' @param size_x (numeric) Font size of the horizontal axis tick labels.
    #' @param size_y (numeric) Font size of the vertical axis tick labels.
    #' @param size_legned (numeric) Font size of the legend.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$pcoa(id_col = 1, group = "group2")
    #' a$pcoa_plot(color_scheme = c("#00b0f6", "#FFC24B", "#f8766d", "#ae876d", "#AFC24B"),
    #'             seed = 123, custom_order = c("A", "B", "R", "D", "S"),
    #'
    #'             size_point = 8, size_point_legend = 8, spacing_legend_point = 1.2,
    #'             spacing_legend_title = 0.5, legend_ncol = 1,
    #'
    #'             title = "PCoA", title_legend = "Group",
    #'
    #'             size_title = 28, size_title_x = 20, size_title_y = 20,
    #'             size_title_legned = 24,
    #'
    #'             size_x = 18, size_y = 18, size_legned = 16, filename = "PCoA",
    #'             file_width = 12, file_height = 9)
    #' }
    #'
    pcoa_plot = function(color_scheme = NULL, custom_order = NULL, seed = 123,
                         size_point = 4.5, size_point_legend = 8, spacing_legend_point = 1.2,
                         spacing_legend_title = 0.5, legend_ncol = 1, label_is = T, size_label = 5,
                         title = "PCoA", title_sub = NULL, title_legend = "Group",
                         size_title = 28, size_title_sub = 16, size_title_x = 20, size_title_y = 20,
                         size_title_legned = 24, size_x = 18, size_y = 18, size_legned = 16,
                         filename = "PCoA", file_width = 12, file_height = 9)
    {
      private$.chart = amplysis::pcoa_plot(
        data = private$.data,
        color_scheme = color_scheme,
        custom_order = custom_order,
        seed = seed,

        size_point = size_point,
        size_point_legend = size_point_legend,
        spacing_legend_point = spacing_legend_point,
        spacing_legend_title = spacing_legend_title,
        legend_ncol = legend_ncol,
        label_is = label_is,
        size_label = size_label,

        title = title,
        title_sub = title_sub,
        title_legend = title_legend,

        size_title = size_title,
        size_title_sub = size_title_sub,
        size_title_x = size_title_x,
        size_title_y = size_title_y,
        size_title_legned = size_title_legned,

        size_x = size_x,
        size_y = size_y,
        size_legned = size_legned,

        filename = filename,
        file_width = file_width,
        file_height = file_height
      )
      print(private$.chart)
    },


    ## NMDS
    ############################################################################
    #' @description
    #' In bioinformatics analysis, Non-Metric Multidimensional Scaling (NMDS) is
    #' used to reduce the dimensionality of high-dimensional data and represent the
    #' similarity or dissimilarity between samples in a lower-dimensional space.
    #' The "nmds()" function utilizes "vegan::vegdist" to calculate the Bray-Curtis
    #' distance and employs "vegan::metaMDS" for conducting NMDS ordination analysis.
    #'
    #' @param id_col (integer) The OTU_ID column is in which column, defaulting to 0
    #' means there is no OTU_ID column, and the data is already numeric.
    #' @param group (Required, character) Grouping information. please enter the
    #' column name of the grouping information in the metadata table.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$nmds(1, "group2")
    #' }
    #'
    nmds = function(id_col = 1, group = "group")
    {
      private$.data = amplysis::nmds(
        otu = private$.otu,
        metadata = private$.metadata,
        id_col = id_col,
        group = group
      )
    },


    ## NMDS_Plot
    #' @description
    #' After processing with the `nmds()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #' @param custom_order Custom legend order.
    #' @param seed Seed
    #'
    #' @param size_point (numeric) The size of points.
    #' @param size_point_legend (numeric) The size of legend points.
    #' @param spacing_legend_point (numeric) The internal padding of the legend.
    #' @param spacing_legend_title (numeric) The spacing between the legend title and the body
    #' @param legend_ncol (integer) Number of columns in the legend.
    #' @param label_is (logical) Whether to display data labels.
    #' @param size_label (numeric) Font size of the data label.
    #'
    #' @param title (character) Main title.
    #' @param title_sub (character) Subtitle.
    #' @param title_legend (character) Legend title.
    #'
    #' @param size_title (numeric) Font size of the main title.
    #' @param size_title_sub (numeric) Font size of the subtitle.
    #' @param size_title_x (numeric) Font size of the horizontal axis title.
    #' @param size_title_y (numeric) Font size of the vertical axis title.
    #' @param size_title_legned (numeric) Font size of legend title.
    #'
    #' @param size_x (numeric) Font size of the horizontal axis tick labels.
    #' @param size_y (numeric) Font size of the vertical axis tick labels.
    #' @param size_legned (numeric) Font size of the legend.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$nmds(1, "group2")
    #' a$nmds_plot(color_scheme = c("#00b0f6", "#FFC24B", "#f8766d", "#ae876d", "#AFC24B"),
    #'             seed = 123, custom_order = c("A", "B", "R", "D", "S"),
    #'
    #'             size_point = 8, size_point_legend = 8, spacing_legend_point = 1.2,
    #'             spacing_legend_title = 0.5, legend_ncol = 1,
    #'
    #'             title = "NMDS", title_legend = "Group",
    #'
    #'             size_title = 28, size_title_x = 20, size_title_y = 20,
    #'             size_title_legned = 24,
    #'
    #'             size_x = 18, size_y = 18, size_legned = 16, filename = "NMDS",
    #'             file_width = 12, file_height = 9)
    #'             }
    #'
    nmds_plot = function(
    color_scheme = NULL,
    custom_order = NULL,
    seed = 123,

    size_point = 4.5,
    size_point_legend = 8,
    spacing_legend_point = 1.2,
    spacing_legend_title = 0.5,
    legend_ncol = 1,
    label_is = T,
    size_label = 5,

    title = "NMDS",
    title_sub = NULL,
    title_legend = "Group",

    size_title = 28,
    size_title_sub = 16,
    size_title_x = 20,
    size_title_y = 20,
    size_title_legned = 24,

    size_x = 18,
    size_y = 18,
    size_legned = 16,

    filename = "NMDS",
    file_width = 12,
    file_height = 9
    )
    {
      private$.chart = amplysis::nmds_plot(
        data = private$.data,
        color_scheme = color_scheme,
        custom_order = custom_order,
        seed = seed,

        size_point = size_point,
        size_point_legend = size_point_legend,
        spacing_legend_point = spacing_legend_point,
        spacing_legend_title = spacing_legend_title,
        legend_ncol = legend_ncol,
        label_is = label_is,
        size_label = size_label,

        title = title,
        title_sub = title_sub,
        title_legend = title_legend,

        size_title = size_title,
        size_title_sub = size_title_sub,
        size_title_x = size_title_x,
        size_title_y = size_title_y,
        size_title_legned = size_title_legned,

        size_x = size_x,
        size_y = size_y,
        size_legned = size_legned,

        filename = filename,
        file_width = file_width,
        file_height = file_height
      )
      print(private$.chart)
    },
    ############################################################################


    ## venn
    #' @description
    #' When the number of samples is less than 5, a Venn diagram can effectively
    #' display the intersections of OTUs among the samples. Although the `venn()`
    #' function allows the `metadata` parameter to be NULL, we still recommend using
    #' the `group` column in the metadata file to control the grouping information.
    #'
    #' It is important to note that when the number of samples is 5 or more, we
    #' recommend using a `upset` diagram to represent the intersections of OTUs.
    #'
    #' @param id_col (integer) The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param group (character) Group 1, please enter the column name of
    #' the grouping information in the metadata table.
    #' @param parallel_method (character) Sample processing methods for the same group:
    #' average, sum, median, none.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$venn(1, "group", "mean")
    #' }
    #'
    venn = function(id_col = 1, group = "group", parallel_method = "mean")
    {
      private$.data = amplysis::venn(
        otu = private$.otu,
        metadata = private$.metadata,
        id_col = id_col,
        group = group,
        parallel_method = parallel_method
      )
    },


    ## venn_plot
    #' @description
    #' After processing with the `venn()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_scheme (character) Color scheme.
    #'
    #' @param show_percentage (logical) Whether to display the percentage of each
    #' set.
    #' @param digits (numeric) Number of decimal places to retain for percentages.
    #'
    #' @param size_set_name (numeric) Title size of each dataset.
    #' @param size_text (numeric) Font size of the text within the dataset.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$venn(1, "group", "mean")
    #' a$venn_plot(color_scheme = c('#aaf200','#0082ff',"#d23aa4","#c777ff", "#79ff79"))
    #' }
    #'
    venn_plot = function(color_scheme = NULL, show_percentage = T, digits = 2,
                         size_set_name = 12, size_text = 7, filename = NULL,
                         file_width = 12, file_height = 9
    )
    {
      private$.chart = amplysis::venn_plot(
        data = private$.data,
        color_scheme = color_scheme,

        show_percentage = show_percentage,
        digits = digits,

        size_set_name = size_set_name,
        size_text = size_text,

        filename = filename,
        file_width = file_width,
        file_height = file_height
      )
      print(private$.chart)
    },
    ############################################################################


    ## Upset
    #' @description
    #' When the number of samples is greater than or equal to 5, the set diagram can
    #' more effectively display the intersection information of OTUs among the
    #' samples. Although the `Upset()` function allows the `metadata` parameter to
    #' be NULL, we still recommend using the `group` column in the metadata file to
    #' control the grouping information.
    #'
    #' @param id_col (integer) The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param group (character) Group 1, please enter the column name of
    #' the grouping information in the metadata table.
    #' @param parallel_method (character) Sample processing methods for the same group:
    #' average, sum, median, none.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$Upset(1, "group", "mean")
    #' }
    #'
    Upset = function(id_col = 1, group = "group", parallel_method = "mean")
    {
      private$.data = amplysis::Upset(
        otu = private$.otu,
        metadata = private$.metadata,
        id_col = id_col,
        group = group,
        parallel_method = parallel_method
      )
    },


    ## Upset_plot
    #' @description
    #' After processing with the `Upset()` function, you can use this function
    #' for visualization.
    #'
    #' @param color_matrix (character) Color scheme of the bar chart in the matrix
    #' plot.
    #' @param n (integer) Number of bars displayed in the bar chart.
    #' @param custom_order (character) Order of the x-axis in the bar chart.
    #' @param order_by (character) Sorting order of the matrix plot, options are
    #' "freq" or "degree".
    #'
    #' @param mb.ratio (character) The proportion of the heights between the bar
    #' plot and the matrix plot.
    #' @param size_point (numeric) The size of the points in the matrix plot.
    #'
    #' @param color_point (character) The color of the points in the matrix plot.
    #' @param color_matrix_shade (character) The color of the shaded areas in the
    #' matrix plot.
    #' @param color_bar (character) The color of the bars in the y-axis bar plot.
    #'
    #' @param title_matrix_x (character) The labels on the x-axis of the matrix plot.
    #' @param title_bar_y (character) The title of the Y-axis in the bar plot.
    #'
    #' @param size_title_matrix_x (numeric) The font size of the x-axis title in the
    #' matrix plot.
    #' @param size_title_matrix_y (numeric) The font size of the y-axis title in the
    #' matrix plot.
    #' @param size_title_bar_y (numeric) The font size of the vertical axis title in
    #' the bar plot.
    #'
    #' @param size_matrix_x (numeric) The size of the tick labels in the matrix plot.
    #' @param size_bar_y (numeric) The font size of the characters on the vertical
    #' axis scale of the bar plot.
    #' @param size_bar_label (numeric) The size of the numbers on the bars in the
    #' bar plot.
    #'
    #' @param queries (a list) Highlight specific elements in the matrix plot or the
    #' bar plot.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @examples
    #' \dontrun{
    #' upset1 <- Upset(otu, metadata, 1, "group", "mean")
    #' Upset_plot(upset1,
    #'            queries = list(
    #'              list(query = intersects,
    #'                   params = list("AA", "AB", "AD"), color="#f06676", active = T),
    #'              list(query = intersects,
    #'                   params = list("AA", "AB"), color="#f06676", active = T)))
    #' }
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$Upset(1, "group", "mean")
    #' a$Upset_plot(
    #'   color_matrix = NULL, n = 40, custom_order = NULL, order_by = "freq",
    #'   mb.ratio = c(0.6, 0.4),
    #'   size_point = 2, color_point = "#505050", color_matrix_shade = "#f89c9f",
    #'   color_bar = '#505050',
    #'   title_matrix_x = "Set Size", title_bar_y = "Intersection Size",
    #'   size_title_matrix_x = 2, size_title_matrix_y = 1.7, size_title_bar_y = 2,
    #'   size_matrix_x = 2, size_bar_y = 2, size_bar_label = 1.75)
    #' }
    #'
    Upset_plot = function(color_matrix = NULL, n = 40, custom_order = NULL,
                          order_by = "freq", mb.ratio = c(0.6, 0.4), size_point = 2,
                          color_point = "#505050", color_matrix_shade = "#f89c9f",
                          color_bar = '#505050', title_matrix_x = "Set Size",
                          title_bar_y = "Intersection Size", size_title_matrix_x = 2,
                          size_title_matrix_y = 1.7, size_title_bar_y = 2, size_matrix_x = 2,
                          size_bar_y = 2, size_bar_label = 1.75, queries = NULL,
                          filename = "Upset", file_width = 16, file_height = 9
    )
    {
      private$.chart = amplysis::Upset_plot(
        data = private$.data,
        color_matrix = color_matrix,

        n = n,
        custom_order = custom_order,
        order_by = order_by,

        mb.ratio = mb.ratio,
        size_point = size_point,

        color_point = color_point,
        color_matrix_shade = color_matrix_shade,
        color_bar = color_bar,

        title_matrix_x = title_matrix_x,
        title_bar_y = title_bar_y,

        size_title_matrix_x = size_title_matrix_x,
        size_title_matrix_y = size_title_matrix_y,
        size_title_bar_y = size_title_bar_y,

        size_matrix_x = size_matrix_x,
        size_bar_y = size_bar_y,
        size_bar_label = size_bar_label,

        queries = queries,

        filename = filename,
        file_width = file_width,
        file_height = file_height
      )
      print(private$.chart)
    },
    ############################################################################


    ## network
    #' @description
    #' The feature table, taxonomy table, and metadata file will be processed into
    #' an edge file and a node file. These two files can be imported into Gephi for
    #' visualization (which we also recommend). Additionally, the network() function
    #' has built-in visualization methods, but they are performance-intensive and
    #' may run slower.
    #'
    #' @param id_col (integer) The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param tax_cla (character) Taxonomic level. Only column names from the tax
    #' table can be entered, for example tax_cla = 'genus'.
    #' @param label (character) The "label" column in the node file generally comes
    #' from a column name in the taxonomy file, with the default being "phylum".
    #' @param group (character) Group 1, please enter the column name of
    #' the grouping information in the metadata table.
    #' @param parallel_method (character) Sample processing methods for the same group:
    #' average, sum, median, none.
    #' @param calc_method (character) Correlation analysis methods available:
    #' "spearman" or "pearson". The default method is "Spearman".
    #' @param cluster_method (integer) Clustering method, options range from 1 to 11, with the
    #' default being 1 (Louvain method).
    #'
    #' A list of all methods:
    #'
    #' "1" = "Louvain (default)",
    #'
    #' "2" = "Edge Betweenness",
    #'
    #' "3" = "Fluid Communities",
    #'
    #' "4" = "Infomap",
    #'
    #' "5" = "Label Propagation",
    #'
    #' "6" = "Leading Eigen",
    #'
    #' "7" = "Leiden",
    #'
    #' "8" = "Optimal",
    #'
    #' "9" = "Spinglass",
    #'
    #' "10" = "Walktrap",
    #'
    #' "11" = "Fast Greedy")
    #'
    #' @param normalize_flag (logical) Standardize the columns "eigenvector", "closeness",
    #' "constraint", and "pagerank" in the node file. Default is TRUE.
    #'
    #' @param .r (numeric) Correlation. Default is 0.6, meaning nodes with an
    #' r-value > 0.6 will be retained.
    #' @param .p (numeric) Significance test p-value. Default is 0.05, meaning nodes with
    #' p-values less than 0.05 will be retained.
    #'
    #' @examples
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$network()}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$network(id_col = 1, .r = 0.6, .p = 0.05)}
    #' \dontrun{
    #' a = amplysis$new(otu, tax, rep, tree, env, metadata)
    #' a$network(id_col = 1, tax_cla = "genus", label = "phylum",
    #'           group = "group", parallel_method = "mean", calc_method = "spearman",
    #'           cluster_method = 1, normalize_flag = TRUE, .r = 0.6, .p = 0.05)}
    #'
    network = function(id_col = 1, tax_cla = "genus", label = "phylum", group = "group",
                       parallel_method = "mean", calc_method = "spearman",
                       cluster_method = 1, normalize_flag = TRUE, .r = 0.6, .p = 0.05)
    {
      private$.data = amplysis::network(
        otu = private$.otu,
        tax = private$.tax,
        metadata = private$.metadata,
        id_col = id_col,
        tax_cla = tax_cla,
        label = label,
        group = group,
        parallel_method = parallel_method,
        calc_method = calc_method,
        cluster_method = cluster_method,
        normalize_flag = normalize_flag,
        .r = .r,
        .p = .p
      )
    }


  ),


  ##
  private = list(),


  ##
  active = list(
    #' @field otu Get or set `OTU/ASV` table
    otu = function(otu) {
      if(missing(otu))
        return(private$.otu)
      private$.otu = otu
    },


    #' @field tax Get or set `tax` table
    tax = function(tax) {
      if(missing(tax))
        return(private$.tax)
      private$.tax = tax
    },


    #' @field rep_seq Get or set `rep_seq`
    rep_seq = function(rep_seq) {
      if(missing(rep_seq))
        return(private$.rep_seq)
      private$.rep_seq = rep_seq
    },


    #' @field tree Get or set `phylogenetic tree` file
    tree = function(tree) {
      if(missing(tree))
        return(private$.tree)
      private$.tree = tree
    },


    #' @field metadata Get or set `metadata` table
    metadata = function(metadata) {
      if(missing(metadata))
        return(private$.metadata)
      private$.metadata = metadata
    },



    #' @field data Get `data`
    data = function() {
      return(private$.data)
    },


    #' @field chart Get `chart`
    chart = function() {
      return(private$.chart)
    }
  )
)


