% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_plot.R
\name{write_plot}
\alias{write_plot}
\title{Save a ggplot or ComplexHeatmap as a standalone file
and return path to the file.
I tries to guess the size of the plot if possible
(downstream of set_panel_size() or after defining a heatmap size)}
\usage{
write_plot(
  x,
  filename,
  device = grDevices::png,
  width = NA,
  height = NA,
  units = "cm",
  ...
)
}
\arguments{
\item{x}{a ggplot or ComplexHeatmap.}

\item{filename}{The basename of the output file}

\item{device}{Device to use. Defaults to PNG.}

\item{width, height, units}{The settings for the output file.}

\item{...}{Arguments passed to device or draw.}
}
\description{
Save a ggplot or ComplexHeatmap as a standalone file
and return path to the file.
I tries to guess the size of the plot if possible
(downstream of set_panel_size() or after defining a heatmap size)
}
