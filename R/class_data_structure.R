#' @title R6 Class: otu_metadata_data_chart
#'
#' @description
#' Data structure:
#' otu,
#' metadata,
#' data,
#' chart.
#'
#' @importFrom R6 R6Class
otu_metadata_data_chart = R6::R6Class(
  # Define an R6 class: Base class.
  "otu_metadata_data_chart",


  public = list(
    #' @description
    #' initialization function.
    #'
    #' @param otu otu table
    #' @param metadata metadata table
    #'
    #' @return R6 class
    #'
    #' @examples \dontrun{ otu_metadata_data_chart$new(otu, metadata) }
    #'
    initialize = function(otu = NA, metadata = NA) {
      private$.otu = otu
      private$.metadata = metadata
      private$.data = NA
      private$.chart = NA
    },


    #' @description
    #' verride the print function.
    #'
    #' @param ...
    #'
    #' @return gg, ggplot2, chart
    #'
    #' @examples
    #' \dontrun{
    #' a <- otu_metadata_data_chart$new(otu, metadata)
    #' print(a)
    #' }
    print = function(...) {
      print(private$.chart)
    }
  ),


  private = list(
    .otu = NA,       # otu table.
    .metadata = NA,  # metadata table.
    .data = NA,      # The data obtained after processing through the data analysis functions.
    .chart = NA      # The chart object obtained after processing through the visualization functions.
  ),


  active = list(
    #' @field otu Get or set `OTU/ASV` table
    otu = function(otu) {
      if(missing(otu))
        return(private$.otu)
      private$.otu = otu
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
################################################################################


#' @title R6 Class: otu_tax_metadata_data_chart
#'
#' @description
#' Data structure:
#' otu,
#' tax,
#' metadata,
#' data,
#' chart.
#'
#' @importFrom R6 R6Class
otu_tax_metadata_data_chart = R6::R6Class(
  # inherit
  # Add a new `tax` attribute.

  "otu_tax_metadata_data_chart",

  #
  inherit = otu_metadata_data_chart,


  ##
  public = list(
    #' @description
    #' initialization function.
    #'
    #' @param otu otu table
    #' @param tax tax table
    #' @param metadata metadata table
    #'
    #' @return R6 class
    #'
    #' @examples \dontrun{ otu_tax_metadata_data_chart$new(otu, tax, metadata) }
    #'
    initialize = function(otu = NA, tax = NA, metadata = NA) {
      private$.otu = otu
      private$.tax = tax
      private$.metadata = metadata
      private$.data = NA
      private$.chart = NA
    }
  ),


  ##
  private = list(
    .tax = NA  # tax table
  ),


  ##
  active = list(
    #' @field tax Get or set `tax` table
    tax = function(tax) {
      if(missing(tax))
        return(private$.tax)
      private$.tax = tax
    }
  )
)
################################################################################


#' @title R6 Class: otu_tax_rep_metadata_data_chart
#'
#' @description
#' Data structure:
#' otu,
#' tax,
#' rep_seq,
#' metadata,
#' data,
#' chart.
#'
#' @importFrom R6 R6Class
otu_tax_rep_metadata_data_chart = R6::R6Class(
  # inherit
  # Add a new `rep_seq` attribute.

  "otu_tax_rep_metadata_data_chart",

  #
  inherit = otu_tax_metadata_data_chart,


  ##
  public = list(
    #' @description
    #' initialization function.
    #'
    #' @param otu otu table
    #' @param tax tax table
    #' @param rep_seq representative sequence
    #' @param metadata metadata table
    #'
    #' @return R6 class
    #'
    #' @examples \dontrun{ otu_tax_rep_metadata_data_chart$new(otu, tax, rep_seq, metadata) }
    #'
    initialize = function(otu = NA, tax = NA, rep_seq = NA, metadata = NA) {
      private$.otu = otu
      private$.tax = tax
      private$.rep_seq = rep_seq
      private$.metadata = metadata
      private$.data = NA
      private$.chart = NA
    }
  ),


  ##
  private = list(
    .rep_seq = NA  # representative sequence
  ),


  ##
  active = list(
    #' @field rep_seq Get or set `rep_seq`
    rep_seq = function(rep_seq) {
      if(missing(rep_seq))
        return(private$.rep_seq)
      private$.rep_seq = rep_seq
    }
  )
)
################################################################################


#' @title R6 Class: otu_tax_rep_tree_metadata_data_chart
#'
#' @description
#' Data structure:
#' otu,
#' tax,
#' rep_seq,
#' tree,
#' metadata,
#' data,
#' chart.
#'
#' @importFrom R6 R6Class
otu_tax_rep_tree_metadata_data_chart = R6::R6Class(
  # inherit
  # Add a new `tree` attribute.

  "otu_tax_rep_tree_metadata_data_chart",

  #
  inherit = otu_tax_rep_metadata_data_chart,


  ##
  public = list(
    #' @description
    #' initialization function.
    #'
    #' @param otu otu table
    #' @param tax tax table
    #' @param rep_seq representative sequence
    #' @param tree phylogenetic tree file
    #' @param metadata metadata table
    #'
    #' @return R6 class
    #'
    #' @examples \dontrun{ otu_tax_rep_tree_metadata_data_chart$new(otu, tax, rep_seq, tree, metadata) }
    #'
    initialize = function(otu = NA, tax = NA, rep_seq = NA, tree = NA, metadata = NA) {
      private$.otu = otu
      private$.tax = tax
      private$.rep_seq = rep_seq
      private$.tree = tree
      private$.metadata = metadata
      private$.data = NA
      private$.chart = NA
    }
  ),


  ##
  private = list(
    .tree = NA  # phylogenetic tree file
  ),


  ##
  active = list(
    #' @field tree Get or set `phylogenetic tree` file
    tree = function(tree) {
      if(missing(tree))
        return(private$.tree)
      private$.tree = tree
    }
  )
)

################################################################################


#' @title R6 Class: otu_tax_rep_tree_env_metadata_data_chart
#'
#' @description
#' Data structure:
#' otu,
#' tax,
#' rep_seq,
#' tree,
#' env,
#' metadata,
#' data,
#' chart.
#'
#' @importFrom R6 R6Class
otu_tax_rep_tree_env_metadata_data_chart = R6::R6Class(
  # inherit
  # Add a new `env` attribute.


  "otu_tax_rep_tree_env_metadata_data_chart",

  #
  inherit = otu_tax_rep_tree_metadata_data_chart,


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
    #' @examples \dontrun{ otu_tax_rep_tree_env_metadata_data_chart$new(otu, tax, rep_seq, tree, env, metadata) }
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
    }
  ),


  ##
  private = list(
    .env = NULL  # environmental factor
  ),


  ##
  active = list(
    #' @field env Get or set `environmental factor` file
    env = function(env) {
      if(missing(env))
        return(private$.env)
      private$.env = env
    }
  )
)


