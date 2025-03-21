% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/km_summary.R
\name{stat_km_summary}
\alias{stat_km_summary}
\title{Calculates the number of observations in each stratum.}
\usage{
stat_km_summary(
  mapping = NULL,
  data = NULL,
  geom = "km",
  position = "identity",
  show.legend = NA,
  inherit.aes = TRUE,
  se = TRUE,
  trans = "identity",
  firstx = 0,
  firsty = 1,
  type = "kaplan-meier",
  start.time = 0,
  range = NULL
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[ggplot2:aes]{aes()}}. If specified and
\code{inherit.aes = TRUE} (the default), it is combined with the default mapping
at the top level of the plot. You must supply \code{mapping} if there is no plot
mapping.}

\item{data}{The data to be displayed in this layer. There are three
options:

If \code{NULL}, the default, the data is inherited from the plot
data as specified in the call to \code{\link[ggplot2:ggplot]{ggplot()}}.

A \code{data.frame}, or other object, will override the plot
data. All objects will be fortified to produce a data frame. See
\code{\link[ggplot2:fortify]{fortify()}} for which variables will be created.

A \code{function} will be called with a single argument,
the plot data. The return value must be a \code{data.frame}, and
will be used as the layer data. A \code{function} can be created
from a \code{formula} (e.g. \code{~ head(.x, 10)}).}

\item{geom}{The geometric object to use to display the data for this layer.
When using a \verb{stat_*()} function to construct a layer, the \code{geom} argument
can be used to override the default coupling between stats and geoms. The
\code{geom} argument accepts the following:
\itemize{
\item A \code{Geom} ggproto subclass, for example \code{GeomPoint}.
\item A string naming the geom. To give the geom as a string, strip the
function name of the \code{geom_} prefix. For example, to use \code{geom_point()},
give the geom as \code{"point"}.
\item For more information and other ways to specify the geom, see the
\link[ggplot2:layer_geoms]{layer geom} documentation.
}}

\item{position}{A position adjustment to use on the data for this layer. This
can be used in various ways, including to prevent overplotting and
improving the display. The \code{position} argument accepts the following:
\itemize{
\item The result of calling a position function, such as \code{position_jitter()}.
This method allows for passing extra arguments to the position.
\item A string naming the position adjustment. To give the position as a
string, strip the function name of the \code{position_} prefix. For example,
to use \code{position_jitter()}, give the position as \code{"jitter"}.
\item For more information and other ways to specify the position, see the
\link[ggplot2:layer_positions]{layer position} documentation.
}}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[ggplot2:borders]{borders()}}.}

\item{firstx, firsty}{the starting point for the survival curves. By default,
the plot program obeys tradition by having the plot start at (0,1).}

\item{...}{Other arguments passed to \code{survival::survfit.formula}}
}
\value{
a data.frame with additional columns: \item{x}{x in data}
  \item{y}{Kaplan-Meier Survival Estimate at x}
  \item{label}{The number of observations as "n = n"}
}
\description{
Calculates the number of observations in each stratum.
}
\section{Aesthetics}{

\code{stat_km_summary} understands the following aesthetics (required aesthetics
are in bold):
\itemize{
  \item \strong{\code{time}} The survival times
  \item \strong{\code{status}} The censoring indicator, see \link[survival]{Surv} for more information.
  \item \code{alpha}
  \item \code{color}
  \item \code{linetype}
  \item \code{size}
}
}

