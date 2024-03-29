% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bubble_heatmap.R
\name{bubble_heatmap}
\alias{bubble_heatmap}
\title{Draw a bubble heatmap.}
\usage{
bubble_heatmap(
  data,
  row,
  column,
  colour,
  size,
  row_key = NULL,
  column_key = NULL,
  dot_size = unit(5, "mm"),
  dot_scale = 0.9,
  contour = TRUE,
  fill = TRUE,
  colour_scale = circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#377EB8",
    "white", "#E41A1C")),
  contour_max_gp = gpar(col = "gray", fill = NA, lty = 1),
  rect_gp = gpar(type = "none"),
  ...
)
}
\arguments{
\item{data}{A tibble containing a "long version" of values and proportions}

\item{row, column}{The variables to be used as columns and rows to represent the matrix}

\item{colour, size}{The variables that will be mapped to the heatmap colour intensity and the dot size}

\item{row_key, column_key}{An optional variable to split the columns and rows (groups)}

\item{dot_size}{The maximum size of the dot/cell of the heatmap.}

\item{dot_scale}{A scaling factor to avoid dots to touch each other (defaults to 0.9).}

\item{contour, fill}{If TRUE uses the colour_scale for the contour and/or the fillin. Otherwise a colour (character).}

\item{colour_scale}{A colour mapping function such as generated by \code{circlize::colorRamp2}.}

\item{contour_max_gp}{Graphical parameter of the maximum contour line (defaults to a gray line)}

\item{rect_gp}{passed to the rect_gp argument of ComplexHeatmap: Should be usualy set to \code{gpar(type = "none")}}

\item{...}{Other arguments passed to \code{ComplexHeatmap::Heatmap}}
}
\value{
a ComplexHeatmap
}
\description{
Draw a bubble heatmap.
}
