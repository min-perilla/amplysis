#' @title Alpha Diversity Analysis
#'
#' @description
#' Alpha diversity analysis includes the calculation of Shannon, Simpson, Chao1,
#' Ace, Pielou, Goods_coverage, and PD indices. It supports significance
#' analysis using one-way ANOVA and multiple comparison methods.
#' The significance marking method uses the letter marking method.
#'
#' @param otu otu table
#' @param metadata metadata table
#' @param id_col (integer) The column number of the OTU ID column in the
#' OTU table, by default, is 1.
#' @param group (character) Group 1, please enter the column name of
#' the grouping information in the metadata table.
#' @param parallel_method (character) Sample processing methods for the same group:
#' average, sum, median, none.
#' @param tree Phylogenetic tree (rooted tree) file.
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
#' @return A list of various alpha diversity indices.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- alpha(otu, metadata, id_col = 1, group = "group",
#' parallel_method = "mean", tree = tree, method = 1)
#' }
#'
#' @importFrom agricolae duncan.test HSD.test LSD.test REGW.test scheffe.test
#' SNK.test waller.test
#' @importFrom dplyr everything group_by group_cols left_join summarise summarize_at
#' @importFrom ggplot2 vars
#' @importFrom picante pd
#' @importFrom stats aov as.formula setNames
#' @importFrom utils write.table
#' @importFrom vegan diversity estimateR specnumber
#'
alpha <- function(otu, metadata, id_col = 1, group = "group",
         parallel_method = "mean", tree = NULL, method = 1)
{
  if (!method %in% 1:7) {
    stop("Invalid method. Please enter a number between 1 and 7.\n",
         "1: Tukey-Hsd\n",
         "2: Fisher-LSD\n",
         "3: S-N-K, Student-Newman-Keuls\n",
         "4: 4 Duncan(new)\n",
         "5: Scheffe\n",
         "6: Waller-Duncan\n",
         "7: REGW\n")
  } else {
    method_name <- switch(as.character(method),
                          "1" = "Turkey-HSD test",
                          "2" = "Fisher-LSD test",
                          "3" = "Student-Newman-Keuls test",
                          "4" = "Duncan test",
                          "5" = "Scheffe test",
                          "6" = "Waller-Duncan test",
                          "7" = "REGW test")
    cat("method: ", method_name, "\n", sep = "")
  }

  #
  if (!all(group %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(
      !all(group %in% colnames(metadata)), "group"),
      "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup: `", group, "`\n\033[30m", sep = "")
  }


  ## 格式检查
  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "parallel" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` column and the `parallel` column!",
         "\nsample: Sample ID (unique)",
         "\nparallel: Parallel sample identifier")
  }


  ##
  metadata2 <- metadata[, c("sample", "parallel", group)]
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    cat("The following row numbers contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No 'NA' values found in the grouping information.")
  }


  ##
  # Obtain the column number of the OTU_ID.
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu <- as.data.frame(otu)            # Convert to data.frame
    rownames(otu) <- otu[, id_col]       # Rename row names
    otu <- otu[, -id_col, drop = FALSE]  # Remove the OTU ID column

  } else {}  # No OTU_ID column


  ##
  allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))
  parallel_method <- base::tolower(parallel_method)

  #
  if(!parallel_method %in% allowedMethods) {
    stop("Please enter the correct argument for the parameter 'parallel_method':\n",
         "Process according to the 'parallel' column in the `metadata` table, ",
         "\nsamples with the same 'parallel' value are considered parallel samples.\n",
         "`mean`  : Calculate the average\n",
         "`sum`   : Calculate the sum\n",
         "`median`: Calculate the median\n",
         "`none`  : Do not process parallel samples\n")
  } else {
    cat("\033[32mParallel parallel_method: `", parallel_method, "`\n\033[30m", sep = "")
  }


  ##
  if(parallel_method %in% c("mean", "sum", "median")) {
    otu_t <- as.data.frame(t(otu))
    otu_t_g <- merge(otu_t, metadata2, by.x = "row.names", by.y = "sample",
                     all.x = TRUE, sort = F)
    row.names(otu_t_g) <- otu_t_g[, 1]
    otu_t_g <- otu_t_g[, -1]
    metadata3 <- otu_t_g[-c(1:ncol(otu_t))]
    rownames(metadata3) <- NULL
    colnames(metadata3)[1] <- "sample"

    #
    otu_t2 <- otu_t %>%
      dplyr::group_by(metadata3[["sample"]]) %>%
      dplyr::summarize_at(ggplot2::vars(-dplyr::group_cols()), parallel_method)

    otu_t2 <- as.data.frame(otu_t2)
    row.names(otu_t2) <- otu_t2[, 1]
    otu_t2 <- otu_t2[, -1]

    #
    otu2 <- as.data.frame(t(otu_t2))
    metadata4 <- unique(as.data.frame(metadata3[["sample"]]))
    colnames(metadata4)[1] <- "sample"
    rownames(metadata4) <- NULL

    matched_indices <- match(metadata4[["sample"]], metadata3[["sample"]])
    unique_group <- metadata3[[group]][matched_indices]
    metadata5 <- data.frame(metadata4, unique_group)
    colnames(metadata5) <- c("sample", "group")

    #
  } else if(parallel_method %in% c("none")) {
    otu2 <- otu
    metadata5 <- metadata2
    colnames(metadata5) <- c("sample", "parallel", "group")
  }


  ##
  Shannon <- as.data.frame(vegan::diversity(otu2, index = "shannon", MARGIN = 2, base = exp(1)))
  Simpson <- as.data.frame(vegan::diversity(otu2, index = "simpson", MARGIN = 2, base =  exp(1)))
  Richness <- as.data.frame(vegan::specnumber(otu2, MARGIN = 2))

  #
  colnames(Shannon)[1] <- "Shannon"
  colnames(Simpson)[1] <- "Simpson"
  colnames(Richness)[1] <- "Richness"

  #
  obs_chao_ace <- as.data.frame(t(vegan::estimateR(t(apply(otu2, 2, floor)))))
  Chao1 <- stats::setNames(as.data.frame(obs_chao_ace[, 2]), "Chao1")
  Ace <- stats::setNames(as.data.frame(obs_chao_ace[, 4]), "Ace")
  obs <- stats::setNames(as.data.frame(obs_chao_ace[, 1]), "obs")

  #
  Pielou <- Shannon / log(Richness)
  colnames(Pielou)[1] <- "Pielou_e"

  #
  Goods_coverage <- as.data.frame(1 - colSums(otu2 == 1) / colSums(otu2))
  colnames(Goods_coverage)[1] <- "Goods_coverage"


  ##
  # PD
  if(!is.null(tree)){
    # ASV_1 -> 'ASV_1'
    otu_pd <- otu2
    row_names <- rownames(otu_pd)
    new_row_names <- paste0("'", row_names, "'")
    rownames(otu_pd) <- new_row_names
    otu_pd_t <- as.data.frame(t(otu_pd))
    PD <- picante::pd(otu_pd_t, tree, include.root = TRUE)
    colnames(PD) <- c("PD", "SR")


    ##
    index <- cbind(Shannon, Simpson, Chao1, Ace, obs, Richness, Pielou,
                   Goods_coverage, PD)



    ###
  } else {
    index <- cbind(Shannon, Simpson, Chao1, Ace, obs, Richness, Pielou,
                   Goods_coverage)
  }


  ##
  sample = rownames(index)
  index <- data.frame(sample, index)
  rownames(index) <- NULL


  ##
  index <- merge(index, metadata5, by = "sample", all.x = TRUE, sort = F)


  ##
  utils::write.table(index, file = "alpha diversity index.csv", row.names = F,
                     col.names = T, sep = ",")


  ##
  Shannon <- index[, c("sample", "Shannon", "group")]
  Simpson <- index[, c("sample", "Simpson", "group")]
  Chao1 <- index[, c("sample", "Chao1", "group")]
  Ace <- index[, c("sample", "Ace", "group")]
  Pielou <- index[, c("sample", "Pielou_e", "group")]
  Goods_coverage <- index[, c("sample", "Goods_coverage", "group")]
  if(!is.null(tree)) {
    PD <- index[, c("sample", "PD", "group")]
  }


  ##
  column_to_rownames <- function(data){
    rownames(data) <- data[, 1]       # Rename row names
    data <- data[, -1, drop = FALSE]  # Remove the OTU ID column\
    return(data)
  }

  Shannon <- column_to_rownames(Shannon)
  Simpson <- column_to_rownames(Simpson)
  Chao1 <- column_to_rownames(Chao1)
  Ace <- column_to_rownames(Ace)
  Pielou <- column_to_rownames(Pielou)
  Goods_coverage <- column_to_rownames(Goods_coverage)
  if(!is.null(tree)) {
    PD <- column_to_rownames(PD)
  }


  ##
  Shannon[["group"]] <- factor(Shannon[["group"]])
  Simpson[["group"]] <- factor(Simpson[["group"]])
  Chao1[["group"]] <- factor(Chao1[["group"]])
  Ace[["group"]] <- factor(Ace[["group"]])
  Pielou[["group"]] <- factor(Pielou[["group"]])
  Goods_coverage[["group"]] <- factor(Goods_coverage[["group"]])
  if(!is.null(tree)) {
    PD[["group"]] <- factor(PD[["group"]])
  }



  ##
  create_formula <- function(data) {
    response_var <- names(data)[1]
    group_var <- names(data)[2]
    formula <- stats::as.formula(paste(response_var, "~", group_var))
    return(formula)
  }

  #
  formula_shannon <- create_formula(Shannon)
  formula_simpson <- create_formula(Simpson)
  formula_chao1 <- create_formula(Chao1)
  formula_ace <- create_formula(Ace)
  formula_pielou <- create_formula(Pielou)
  formula_goods <- create_formula(Goods_coverage)
  if(!is.null(tree)) {
    formula_pd <- create_formula(PD)
  }


  #
  Shannon_aov <- stats::aov(formula_shannon, data = Shannon)
  Simpson_aov <- stats::aov(formula_simpson, data = Simpson)
  Chao1_aov <- stats::aov(formula_chao1, data = Chao1)
  Ace_aov <- stats::aov(formula_ace, data = Ace)
  Pielou_aov <- stats::aov(formula_pielou, data = Pielou)
  Goods_coverage_aov <- stats::aov(formula_goods, data = Goods_coverage)
  if(!is.null(tree)) {
    PD_aov <- stats::aov(formula_pd, data = PD)
  }


  #
  print(summary(Shannon_aov))
  cat("Pr(>F) value below the 0.05 level indicates significant differences
      (overall level)\n\n")
  print(summary(Simpson_aov))
  cat("Pr(>F) value below the 0.05 level indicates significant differences
      (overall level)\n\n")
  print(summary(Chao1_aov))
  cat("Pr(>F) value below the 0.05 level indicates significant differences
      (overall level)\n\n")
  print(summary(Ace_aov))
  cat("Pr(>F) value below the 0.05 level indicates significant differences
      (overall level)\n\n")
  print(summary(Pielou_aov))
  cat("Pr(>F) value below the 0.05 level indicates significant differences
      (overall level)\n\n")
  print(summary(Goods_coverage_aov))
  cat("Pr(>F) value below the 0.05 level indicates significant differences
      (overall level)\n\n")
  if(!is.null(tree)) {
    print(summary(PD_aov))
    cat("Pr(>F) value below the 0.05 level indicates significant differences
        (overall level)\n\n")
  }


  ##
  multiCompare <- function(data, method) {
    result <- NULL
    sum(is.na(data))
    switch(as.character(method),
           "1" = {
             # Tukey-HSD
             hsd <- agricolae::HSD.test(y = data, trt = "group", alpha = 0.05)
             print(hsd$groups)
             result <- hsd$groups
           },
           "2" = {
             # Fisher-LSD
             lsd <- agricolae::LSD.test(y = data, trt = "group", alpha = 0.05)
             print(lsd$groups)
             result <- lsd$groups
           },
           "3" = {
             # S-N-K, Student-Newman-Keuls
             snk <- agricolae::SNK.test(y = data, trt = "group", alpha = 0.05)
             print(snk$groups)
             result <- snk$groups
           },
           "4" = {
             # Duncan (new)
             dc <- agricolae::duncan.test(y = data, trt = "group", alpha = 0.05)
             print(dc$groups)
             result <- dc$groups
           },
           "5" = {
             # Scheffe
             scheffe <- agricolae::scheffe.test(y = data, trt = "group", alpha = 0.05)
             print(scheffe$groups)
             result <- scheffe$groups
           },
           "6" = {
             # Waller-Duncan
             waller <- agricolae::waller.test(y = data, trt = "group")
             print(waller$groups)
             result <- waller$groups
           },
           "7" = {
             # REGW
             regw <- agricolae::REGW.test(y = data, trt = "group", alpha = 0.05)
             print(regw$groups)
             result <- regw$groups
           },
           {
             stop("Invalid method. Please enter a number between 1 and 7.")
           }
    )

    return(result)
  }

  #
  info_Shannon <- multiCompare(Shannon_aov, 1)
  info_Simpson <- multiCompare(Simpson_aov, 1)
  info_Chao1 <- multiCompare(Chao1_aov, 1)
  info_Ace <- multiCompare(Ace_aov, 1)
  info_Pielou <- multiCompare(Pielou_aov, 1)
  info_Goods_coverage <- multiCompare(Goods_coverage_aov, 1)
  if(!is.null(tree)) {
    info_PD <- multiCompare(PD_aov, 1)
  }



  ##
  get_differ <- function(info)
  {
    group_info <- rownames(info)
    new_data <- cbind(group_info, info)
    new_data <- new_data[, -2]
    colnames(new_data) <- c("group", "differ")
    rownames(new_data) <- NULL

    return(new_data)
  }

  #
  differ_Shannon <- get_differ(info_Shannon)
  differ_Simpson <- get_differ(info_Simpson)
  differ_Chao1 <- get_differ(info_Chao1)
  differ_Ace <- get_differ(info_Ace)
  differ_Pielou <- get_differ(info_Pielou)
  differ_Goods_coverage <- get_differ(info_Goods_coverage)
  if(!is.null(tree)) {
    differ_PD <- get_differ(info_PD)
  }


  ##
  differ_position <- function(index_table, differ, type)
  {
    ##
    index_table2 <- index_table
    index_table2$"sample" <- rownames(index_table2)
    rownames(index_table2) <- NULL
    index_table2 <- index_table2 %>% select("sample", dplyr::everything())
    merged_table <- dplyr::left_join(index_table2, differ, by = c("group" = "group"))


    ##
    max = max(merged_table[, 2])
    min = min(merged_table[, 2])
    x = merged_table[, c("group", type)]
    group_x <- "group"
    differ_y = x %>%
      dplyr::group_by(dplyr::select(x, tidyr::all_of(group_x))) %>%
      dplyr::summarise(max_value = max(get(type)))

    differ_y = as.data.frame(differ_y)
    rownames(differ_y) = differ_y$group
    colnames(differ_y) <- c("group", type)
    merged_table$differ_y <- differ_y[as.character(merged_table$group), ][[type]] + (max - min) * 0.1
    rownames(merged_table) <- merged_table[, 1]
    merged_table <- merged_table[, -1]


    ##
    return(merged_table)
  }

  #
  df_Shannon <- differ_position(Shannon, differ_Shannon, "Shannon")
  df_Simpson <- differ_position(Simpson, differ_Simpson, "Simpson")
  df_Chao1 <- differ_position(Chao1, differ_Chao1, "Chao1")
  df_Ace <- differ_position(Ace, differ_Ace, "Ace")
  df_Pielou <- differ_position(Pielou, differ_Pielou, "Pielou_e")
  df_Goods_coverage <- differ_position(Goods_coverage, differ_Goods_coverage, "Goods_coverage")
  if(!is.null(tree)) {
    df_PD <- differ_position(PD, differ_PD, "PD")
  }



  ##
  result <- NULL

  if(is.null(tree)) {
    result <- list(
      shannon = df_Shannon,
      simpson = df_Simpson,
      chao1 = df_Chao1,
      ace = df_Ace,
      pielou = df_Pielou,
      goods_coverage = df_Goods_coverage)
  } else {
    result <- list(
      shannon = df_Shannon,
      simpson = df_Simpson,
      chao1 = df_Chao1,
      ace = df_Ace,
      pielou = df_Pielou,
      goods_coverage = df_Goods_coverage,
      pd = df_PD)
  }

  return(result)
}
