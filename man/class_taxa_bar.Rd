% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_taxa_bar.R
\name{class_taxa_bar}
\alias{class_taxa_bar}
\title{R6 Class: class_taxa_bar}
\value{
class_taxa_bar: otu, tax, metadata, data, chart
}
\description{
The class_taxa_bar class encapsulates the species composition analysis
function \code{taxa_bar()} and its accompanying visualization function
\code{taxa_bar_plot()}, enabling convenient analysis.
}
\examples{

## ------------------------------------------------
## Method `class_taxa_bar$plot`
## ------------------------------------------------

\dontrun{
a <- class_taxa_bar$new(otu, tax, metadata)
a$plot()
}

}
\section{Super classes}{
\code{\link[amplysis:otu_metadata_data_chart]{amplysis::otu_metadata_data_chart}} -> \code{\link[amplysis:otu_tax_metadata_data_chart]{amplysis::otu_tax_metadata_data_chart}} -> \code{class_taxa_bar}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-class_taxa_bar-plot}{\code{class_taxa_bar$plot()}}
\item \href{#method-class_taxa_bar-clone}{\code{class_taxa_bar$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="amplysis" data-topic="otu_metadata_data_chart" data-id="print"><a href='../../amplysis/html/otu_metadata_data_chart.html#method-otu_metadata_data_chart-print'><code>amplysis::otu_metadata_data_chart$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="amplysis" data-topic="otu_tax_metadata_data_chart" data-id="initialize"><a href='../../amplysis/html/otu_tax_metadata_data_chart.html#method-otu_tax_metadata_data_chart-initialize'><code>amplysis::otu_tax_metadata_data_chart$initialize()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-class_taxa_bar-plot"></a>}}
\if{latex}{\out{\hypertarget{method-class_taxa_bar-plot}{}}}
\subsection{Method \code{plot()}}{
The \code{plot()} function integrates the \code{taxa_bar()} and \code{taxa_bar_plot()}
functions, allowing for one-click plotting of species stacked plots.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{class_taxa_bar$plot(
  id_col = 1,
  tax_cla = "phylum",
  group1 = "group",
  group2 = NULL,
  parallel_method = "mean",
  row_n = 8,
  color_scheme = NULL,
  y_abun = "abun",
  custom_order = NULL,
  custom_order_F = NULL,
  bar_type = "fill",
  bar_width = 0.7,
  grid_line = F,
  size_point_legend = 0.75,
  spacing_legend_point = 0.75,
  spacing_legend_title = 0.5,
  legend_ncol = 1,
  title = NULL,
  title_x = "Groups",
  title_y = NULL,
  title_legend = NULL,
  size_title = 28,
  size_title_x = 20,
  size_title_y = 20,
  size_title_legned = 24,
  size_title_facet = 32,
  size_x = 14,
  size_y = 18,
  size_legned = 16,
  filename = NULL,
  file_width = 16,
  file_height = 9
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id_col}}{(integer)The column number of the OTU ID column in the
OTU table, by default, is 1.}

\item{\code{tax_cla}}{(character) Taxonomic level. Only column names from the tax
table can be entered, for example tax_cla = 'genus'.}

\item{\code{group1}}{(Required, character) Group 1, please enter the column name of
the grouping information in the metadata table.}

\item{\code{group2}}{(Optional, character) Group 2 for facetting plots, please enter
the column name of the grouping information in the metadata table.}

\item{\code{parallel_method}}{(character) Parallel sample processing method,
defaulting to mean. Options: mean (average), sum (summation), median (median).}

\item{\code{row_n}}{(integer) Preserve the top N taxa (including the Nth) based on
abundance, while merging taxa with lower abundance into "others".}

\item{\code{color_scheme}}{(character) Color scheme. Please enter a hexadecimal
format color scheme character vector, such as color_scheme = c("#7FC97F",
"#BEAED4", "#FDC086", "#FFFF99", "#386CB0"). If NULL, the default color
scheme will be used.}

\item{\code{y_abun}}{(character) Abundance information. Please enter the column name
for "Abundance Information" in metadata.}

\item{\code{custom_order}}{(character) Customize the order of horizontal axis grouping.}

\item{\code{custom_order_F}}{(character) (Facet plot) Customize the order of
horizontal axis grouping.}

\item{\code{bar_type}}{(character) Type of species stacked plot. Options: "fill"
(relative abundance) or "stack" (absolute abundance).}

\item{\code{bar_width}}{(numeric) Width of each column in the stacked plot.}

\item{\code{grid_line}}{(logical) Grid lines on the background of the stacked plot.}

\item{\code{size_point_legend}}{(numeric) Size of legend points.}

\item{\code{spacing_legend_point}}{(numeric) The spacing inside the legend.}

\item{\code{spacing_legend_title}}{(numeric) The spacing between the legend title
and the main text.}

\item{\code{legend_ncol}}{(integer) Number of columns in the legend.}

\item{\code{title}}{(character) Main title.}

\item{\code{title_x}}{(character) Horizontal axis title.}

\item{\code{title_y}}{(character) Vertical axis title.}

\item{\code{title_legend}}{(character) Legend title.}

\item{\code{size_title}}{(numeric) Font size of the main title.}

\item{\code{size_title_x}}{(numeric) Font size of the horizontal axis title.}

\item{\code{size_title_y}}{(numeric) Font size of the vertical axis title.}

\item{\code{size_title_legned}}{(numeric) Font size of legend title.}

\item{\code{size_title_facet}}{(numeric) Font size of the facet plot titles.}

\item{\code{size_x}}{(numeric) Font size of the horizontal axis tick labels.}

\item{\code{size_y}}{(numeric) Font size of the vertical axis tick labels.}

\item{\code{size_legned}}{(numeric) Font size of the legend.}

\item{\code{filename}}{(character) File name for saving.}

\item{\code{file_width}}{(numeric) Width of the image.}

\item{\code{file_height}}{(numeric) Height of the image.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
data and chart.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
a <- class_taxa_bar$new(otu, tax, metadata)
a$plot()
}

}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-class_taxa_bar-clone"></a>}}
\if{latex}{\out{\hypertarget{method-class_taxa_bar-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{class_taxa_bar$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
