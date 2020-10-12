NULL

#' @importFrom ggplot2 aes ggproto Stat
#' @import ggplot2
#' @importFrom survival Surv survfit.formula
#' @importFrom grid pointsGrob nullGrob unit gpar gList
#' @rdname stat_km
#' @export
StatKmSummary <- ggplot2::ggproto("StatKmSummary", ggplot2::Stat,

                                  compute_group = function(data, scales, firstx = 0, firsty = 1,
                                                           type = "kaplan-meier", start.time = 0, range = NULL, y_value = NULL) {

                                    sf <- survival::survfit.formula(survival::Surv(data$time, data$status) ~ 1, se.fit = FALSE,
                                                                    type = type, start.time = start.time)

                                    if(is.null(sf$surv)) {
                                      x <- rep(sf$time, 2)
                                      sf$surv <- rep(1, length(x))
                                    }

                                    if(is.null(range)) {
                                      my_range <- length(sf$time)
                                    } else {
                                      my_range <- range %in% sf$time
                                    }

                                    if (!is.null(y_value)) sf$surv[my_range] <- max(sf$surv[my_range])

                                    df.out <- data.frame(time = sf$time[my_range], survival = sf$surv[my_range], n = nrow(data))
                                    #df.out <- data.frame(time = sf$time, survival = sf$surv, n = nrow(data))
                                    df.out
                                  },

                                  default_aes = ggplot2::aes(y = ..survival.., x = ..time.., label = paste("n =", ..n..)),
                                  required_aes = c("time", "status")
)

#' Calculates the number of observations in each stratum.
#'
#' @section Aesthetics:
#' \code{stat_km_summary} understands the following aesthetics (required aesthetics
#' are in bold):
#' \itemize{
#'   \item \strong{\code{time}} The survival times
#'   \item \strong{\code{status}} The censoring indicator, see \link[survival]{Surv} for more information.
#'   \item \code{alpha}
#'   \item \code{color}
#'   \item \code{linetype}
#'   \item \code{size}
#' }
#'
#' @inheritParams ggplot2::stat_identity
#' @param firstx,firsty the starting point for the survival curves. By default,
#'   the plot program obeys tradition by having the plot start at (0,1).
#' @param ... Other arguments passed to \code{survival::survfit.formula}
#' @return a data.frame with additional columns: \item{x}{x in data}
#'   \item{y}{Kaplan-Meier Survival Estimate at x}
#'   \item{label}{The number of observations as "n = n"}
#' @export
stat_km_summary <- function(mapping = NULL, data = NULL, geom = "km",
                            position = "identity", show.legend = NA, inherit.aes = TRUE,
                            se = TRUE, trans = "identity", firstx = 0, firsty = 1,
                            type = "kaplan-meier", start.time = 0, range = NULL) {
  ggplot2::layer(
    stat = StatKmSummary,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(trans = trans, firstx = firstx, firsty = firsty,
                  type = type, start.time = start.time, range = range)
  )

}