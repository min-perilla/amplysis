#' @title Co-occurrence Network Analysis
#'
#' @description
#' The feature table, taxonomy table, and metadata file will be processed into
#' an edge file and a node file. These two files can be imported into Gephi for
#' visualization (which we also recommend). Additionally, the network() function
#' has built-in visualization methods, but they are performance-intensive and
#' may run slower.
#'
#' @param otu otu table
#' @param tax taxonomy table
#' @param metadata metadata file
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
#' @return Two CSV files: "edge.csv" and "node.csv".
#' @export
#'
#' @examples
#' \dontrun{
#' network(otu, tax, metadata)}
#' \dontrun{
#' network(otu, tax, metadata, id_col = 1, .r = 0.6, .p = 0.05)}
#' \dontrun{
#' network(otu, tax, metadata, id_col = 1, tax_cla = "genus", label = "phylum",
#'         group = "group", parallel_method = "mean", calc_method = "spearman",
#'         cluster_method = 1, normalize_flag = TRUE, .r = 0.6, .p = 0.05)}
#'
#' @importFrom dplyr across arrange desc filter group_by group_cols mutate
#' rename select ungroup vars where
#'
#' @importFrom stats cor cor.test na.omit sd
#'
#' @importFrom utils install.packages write.csv
#'
#' @importFrom igraph as_adjacency_matrix as_data_frame betweenness closeness
#' constraint cluster_edge_betweenness cluster_fast_greedy
#' cluster_fluid_communities cluster_infomap cluster_label_prop
#' cluster_leading_eigen cluster_leiden cluster_louvain cluster_optimal
#' cluster_spinglass cluster_walktrap degree eigen_centrality
#' graph_from_data_frame membership modularity neighbors page_rank transitivity
#' V vcount vertex_attr
#'
#'
network = function(otu, tax, metadata, id_col = 1, tax_cla = "genus",
                   label = "phylum", group = "group", parallel_method = "mean",
                   calc_method = "spearman", cluster_method = 1, normalize_flag = TRUE,
                   .r = 0.6, .p = 0.05)
{
  valid_methods <- c("spearman", "pearson")

  if (is.na(cluster_method) || cluster_method < 1 || cluster_method > 11) {
    cat("Invalid input. Default method (Louvain) will be used.\n")
    cluster_method <- 1
  }

  if (!(calc_method %in% valid_methods)) {
    stop("Error: Calculation method calc_method must be either 'spearman' or 'pearson'.")
  } else {
    cat("\033[32m", "calc_method: ", calc_method, "\033[39m\n", sep = "")
  }


  ##
  otu_metadata = amplysis::parallel(
    otu = otu, metadata = metadata, id_col = id_col, group = group,
    parallel_method = parallel_method, digits = 0, metadata_out = T)

  otu2 = otu_metadata[["otu"]]
  metadata2 = otu_metadata[["metadata"]]

  otu_tax = amplysis::align_otu_tax(tax, otu2, by.x = 1, by.y = "row.names")
  tax2 = otu_tax[["tax"]]

  row.names(tax2) = tax2[, 1]
  tax2 = tax2[, -1]
  tax2[is.na(tax2)] <- "unknown"


  ##
  otu3 <- otu2 %>%
    dplyr::group_by(tax2[[tax_cla]]) %>%
    dplyr::summarize_at(dplyr::vars(-dplyr::group_cols()), sum) %>%
    dplyr::arrange(dplyr::desc(rowSums(dplyr::across(dplyr::where(is.numeric))))) %>%
    dplyr::ungroup()

  otu3 = as.data.frame(otu3)
  row.names(otu3) = otu3[, 1]
  otu3 = otu3[, -1]

  zero_sum_rows <- rowSums(otu3 == 0) == ncol(otu3)
  otu3_clean <- otu3[!zero_sum_rows, ]

  otu4 <- as.data.frame(t(otu3_clean))


  ##############################################################################
  correlate2 <- function(x, method = "pearson", use = "everything", ...)
  {
    if (!is.matrix(x))
      x <- as.matrix(x)
    y = x
    m = n = ncol(x)
    r <- stats::cor(x, y, use = use, method = method)
    p <- matrix(NA, ncol = m, nrow = n, dimnames = list(rownames(r), colnames(r)))
    id <- expand.grid(1:n, 1:m)
    id <- id[id$Var1 > id$Var2, , drop = FALSE]

    if (!requireNamespace("purrr", quietly = TRUE)) {
      utils::install.packages("purrr")
    }


    purrr::walk2(id$Var1, id$Var2, function(.idx, .idy) {
      tmp <- suppressWarnings(stats::cor.test(x = x[, .idx],
                                              y = y[, .idy],
                                              method = method,
                                              ...))
      p[c(.idx, .idy), c(.idy, .idx)] <<- tmp$p.value
    })
    diag(p) <- 0
    structure(
      .Data = list(
        r = r,
        p = p),
      class = "correlate"
    )
  }


  ##############################################################################
  as_upper_tri <- function(cor_result) {
    r <- cor_result$r
    p <- cor_result$p
    upper_triangle_r <- r
    upper_triangle_p <- p
    upper_triangle_r[lower.tri(r, diag = TRUE)] <- NA
    upper_triangle_p[lower.tri(p, diag = TRUE)] <- NA

    upper_triangle_r <- as.data.frame(as.table(upper_triangle_r))
    upper_triangle_p <- as.data.frame(as.table(upper_triangle_p))

    upper_triangle <- merge(upper_triangle_r, upper_triangle_p, by = c("Var1", "Var2"), sort = F)
    colnames(upper_triangle) <- c(".rownames", ".colnames", "r", "p")

    upper_triangle <- na.omit(upper_triangle)
    row.names(upper_triangle) = NULL

    return(upper_triangle)
  }


  ##############################################################################
  # error1

  rename_columns <- function(data, old_name, new_name) {
    if (old_name %in% names(data)) {
      data <- data %>%
        rename(!!new_name := !!sym(old_name))
    } else {
      stop(paste("The column name", old_name, "does not exist"))
    }
    return(data)
  }

  row_col_name <- ".rownames"
  col_col_name <- ".colnames"

  r <- sym("r")
  p <- sym("p")

  edge <- otu4 %>%
    correlate2(method = "spearman") %>%
    as_upper_tri() %>%
    dplyr::filter(abs(!!r) > .r, !!p < .p) %>%

    dplyr::mutate(type = "Undirected",
                  id = seq_len(dplyr::n()),
                  label = "",
                  sign = ifelse(!!r > 0, "P", "N"),
                  abs_r = abs(!!r)) %>%
    rename_columns(row_col_name, "source") %>%
    rename_columns(col_col_name, "target") %>%
    dplyr::select("source", "target", "type", "id", "label", "r", "p", "sign", "abs_r")


  ##
  {
    edge_igraph <- igraph::graph_from_data_frame(edge, directed = FALSE)
    node <- igraph::as_data_frame(edge_igraph, "vertices")
    }


  ##
  tax3 = tax2[, c(tax_cla, label)]
  tax3 = tax3[!duplicated(tax3[[tax_cla]]), ]
  node = merge(node, tax3, by.x = "name", by.y = "genus", sort = F)
  colnames(node)[2] = "label"


  ##
  {
    Degree <- as.data.frame(igraph::degree(edge_igraph, mode = "all", loops = TRUE))
    Degree <- data.frame(row.names(Degree), Degree, row.names = NULL)
    colnames(Degree) = c("name", "degree")
    node <- merge(node, Degree, by = c("name" = "name"), sort = F)
  }

  {
    Betweenness <- as.data.frame(igraph::betweenness(edge_igraph))
    Betweenness <- data.frame(row.names(Betweenness), Betweenness, row.names = NULL)
    colnames(Betweenness) = c("name", "betweenness")
    node <- merge(node, Betweenness, by = c("name" = "name"), sort = F)
  }

  {
    Eigenvector <- as.data.frame(igraph::eigen_centrality(edge_igraph)$vector)
    Eigenvector <- data.frame(row.names(Eigenvector), Eigenvector, row.names = NULL)
    colnames(Eigenvector) = c("name", "eigenvector")
    node <- merge(node, Eigenvector, by = c("name" = "name"), sort = F)
  }

  {
    Closeness <- as.data.frame(igraph::closeness(edge_igraph))
    Closeness <- data.frame(row.names(Closeness), Closeness, row.names = NULL)
    colnames(Closeness) = c("name", "closeness")
    node <- merge(node, Closeness, by = c("name" = "name"), sort = F)
  }

  {
    Constraint <- as.data.frame(igraph::constraint(edge_igraph))
    Constraint <- data.frame(row.names(Constraint), Constraint, row.names = NULL)
    colnames(Constraint) = c("name", "constraint")
    node <- merge(node, Constraint, by = c("name" = "name"), sort = F)
  }

  {
    PageRank <- as.data.frame(igraph::page_rank(edge_igraph)$vector)
    PageRank <- data.frame(row.names(PageRank), PageRank, row.names = NULL)
    colnames(PageRank) = c("name", "pagerank")
    node <- merge(node, PageRank, by = c("name" = "name"), sort = F)
  }

  {
    ClusteringCoefficient <- as.data.frame(igraph::transitivity(edge_igraph, type = "local"))
    ClusteringCoefficient <- data.frame(row.names(ClusteringCoefficient), ClusteringCoefficient, row.names = NULL)
    colnames(ClusteringCoefficient) = c("name", "clustering_coefficient")
    node <- merge(node, ClusteringCoefficient, by = c("name" = "name"), sort = F)
  }

  cat("You can choose from up to 11 cluster detection methods by entering the corresponding number:\n",
      "1: Louvain (default)\n",
      "2: Edge Betweenness\n",
      "3: Fluid Communities\n",
      "4: Infomap\n",
      "5: Label Propagation\n",
      "6: Leading Eigen\n",
      "7: Leiden\n",
      "8: Optimal\n",
      "9: Spinglass\n",
      "10: Walktrap\n",
      "11: Fast Greedy\n")

  cat("\033[32m", "You are currently using method: ", cluster_method,
      " - ", switch(cluster_method,
                    "1" = "Louvain (default)",
                    "2" = "Edge Betweenness",
                    "3" = "Fluid Communities",
                    "4" = "Infomap",
                    "5" = "Label Propagation",
                    "6" = "Leading Eigen",
                    "7" = "Leiden",
                    "8" = "Optimal",
                    "9" = "Spinglass",
                    "10" = "Walktrap",
                    "11" = "Fast Greedy"),
      "\033[39m\n", sep = "")

  Community = switch(
    cluster_method,
    "1" = igraph::cluster_louvain(edge_igraph),
    "2" = igraph::cluster_edge_betweenness(edge_igraph),
    "3" = igraph::cluster_fast_greedy(edge_igraph),
    "4" = igraph::cluster_fluid_communities(edge_igraph, 2),
    "5" = igraph::cluster_infomap(edge_igraph),
    "6" = igraph::cluster_label_prop(edge_igraph),
    "7" = igraph::cluster_leading_eigen(edge_igraph),
    "8" = igraph::cluster_leiden(edge_igraph),
    "9" = igraph::cluster_optimal(edge_igraph),
    "10" = igraph::cluster_spinglass(edge_igraph),
    "11" = igraph::cluster_walktrap(edge_igraph)
  )

  Community2 <- igraph::membership(Community)
  Community3 <- as.data.frame(Community2)
  Community3 <- data.frame(row.names(Community3), Community3, row.names = NULL)
  colnames(Community3) = c("name", "community")
  node <- merge(node, Community3, by = c("name" = "name"), sort = F)

  Modularity = igraph::modularity(Community)
  cat("The modularity of the network using Louvain is:", Modularity, "\n")


  ##
  if(isTRUE(normalize_flag)){
    normalize_by_column_names <- function(data, columns_to_normalize) {
      for (col in columns_to_normalize) {
        if (col %in% colnames(data)) {
          data[[col]] <- (data[[col]] - min(data[[col]])) / (max(data[[col]]) - min(data[[col]])) * 100
        } else {
          cat("Column", col, "not found in the dataset.\n")
        }
      }
      return(data)
    }

    columns_to_normalize <- c("eigenvector", "closeness", "constraint", "pagerank")
    node <- normalize_by_column_names(node, columns_to_normalize)
  }


  ##
  edge_igraph2 = edge_igraph
  z = Ki = rep.int(0, igraph::vcount(edge_igraph2))
  names(z) = names(Ki) = igraph::V(edge_igraph2)$name

  for (i in V(edge_igraph2)) {
    node_id <- as.numeric(i)
    comm_i <- Community2[node_id]
    nodes_in_comm <- which(Community2 == comm_i)
    Ki[node_id] = length(
      intersect(igraph::neighbors(edge_igraph2, node_id), nodes_in_comm))
  }


  ##
  N <- max(Community2)
  nS <- tabulate(Community2)

  Ksi <- rep(0, max(Community2))
  sigKsi <- rep(0, max(Community2))
  S = NULL


  ##
  for (S in seq_len(max(Community2))) {
    x <- Ki[Community2 == S]

    if (length(x) > 0) {
      Ksi[S] = mean(unlist(x))
      sigKsi[S] = stats::sd(unlist(x))
    } else {
      Ksi[S] = 0
      sigKsi[S] = 0
    }
  }



  ###
  z = (Ki - Ksi[Community2]) / sigKsi[Community2]


  ##
  z[is.infinite(z)] = 0
  z[is.nan(z)] = 0
  Zi = z

  Zi <- data.frame(Ki, Zi, row.names = names(Ki))



  ###
  igraph::V(edge_igraph2)$module = Community2
  memb <- igraph::vertex_attr(edge_igraph2, "module")
  A <- as.data.frame(igraph::as_adjacency_matrix(edge_igraph2, sparse = FALSE))
  Ki2 = colSums(A)
  Ki_sum = t(rowsum(A, memb))

  Pi = 1 - ((1 / Ki2^2) * rowSums(Ki_sum^2))
  Pi = as.data.frame(Pi)


  ##
  Zi_Pi = merge(Zi, Pi, by = "row.names", sort = F)
  Zi_Pi <- na.omit(Zi_Pi)


  ##
  if(!nrow(Zi_Pi) == 0){
    Zi_Pi[which(Zi_Pi$Zi < 2.5 & Zi_Pi$Pi < 0.62),'type'] <- 'Peripherals'
    Zi_Pi[which(Zi_Pi$Zi < 2.5 & Zi_Pi$Pi >= 0.62),'type'] <- 'Connectors'
    Zi_Pi[which(Zi_Pi$Zi >= 2.5 & Zi_Pi$Pi < 0.62),'type'] <- 'Module hubs'
    Zi_Pi[which(Zi_Pi$Zi >= 2.5 & Zi_Pi$Pi >= 0.62),'type'] <- 'Network hubs'
  }

  node2 = merge(node, Zi_Pi, by.x = "name", by.y = "Row.names", sort = F)


  ##
  utils::write.csv(x = edge, file = "edge.csv", row.names = F)
  names(node2)[names(node2) == "name"] <- "ID"
  utils::write.csv(x = node2, file = "node.csv",  row.names = F)

  cat("\033[32m--- Please use the `network_plot()` function for visualization. ---\n\033[0m")

  ##
  result = NULL
  result = list(
    edge = edge,
    node = node2
  )
  return(result)
}
