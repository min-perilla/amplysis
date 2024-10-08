% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/venn_plot.R
\name{venn_plot}
\alias{venn_plot}
\title{Visualization of Venn}
\usage{
venn_plot(
  data,
  color_scheme = NULL,
  show_percentage = T,
  digits = 2,
  size_set_name = 12,
  size_text = 7,
  filename = NULL,
  file_width = 12,
  file_height = 9
)
}
\arguments{
\item{data}{Plotting data.}

\item{color_scheme}{(character) Color scheme.}

\item{show_percentage}{(logical) Whether to display the percentage of each
set.}

\item{digits}{(numeric) Number of decimal places to retain for percentages.}

\item{size_set_name}{(numeric) Title size of each dataset.}

\item{size_text}{(numeric) Font size of the text within the dataset.}

\item{filename}{(character) File name for saving.}

\item{file_width}{(numeric) Width of the image.}

\item{file_height}{(numeric) Height of the image.}
}
\value{
ggvenn plot
}
\description{
After processing with the \code{venn()} function, you can use this function
for visualization.
}
\examples{
\dontrun{
venn1 <- venn(otu, metadata, 1, "group")
venn_plot(venn1, color_scheme = c('#aaf200','#0082ff',"#d23aa4","#c777ff", "#79ff79"))
}

}
