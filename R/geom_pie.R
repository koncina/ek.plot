NULL

#' @export
draw_key_pie <- function(data, params, size) {
  pieGrob(x = 0.5, y = 0.5,
          r = 0.4, prop = 0.25,
          pie_gp = gpar(
            col = NA,
            fill = alpha(data$colour, data$alpha)),
          background_gp = gpar(
            col = NA,
            fill = alpha(data$fill, data$alpha)
          )
  )
}





#' Draw a piechart for binary proportions.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom tidyselect eval_select
#' @importFrom rlang enquo as_name
#' @importFrom ComplexHeatmap Heatmap
#'
#' @export
GeomPie <- ggproto("GeomPie", Geom,

                   required_aes = c("x", "r",  "prop"),

                   non_missing_aes = c("fill", "colour", "y", "npcy"),

                   default_aes = aes(
                     colour = "black",
                     fill = "gray",
                     alpha = 1,
                     npcy = 0.8
                   ),

                   draw_key = draw_key_pie,

                   draw_panel = function(data, panel_params, coord) {
                     test <<- data
                     coords <- coord$transform(data, panel_params)

                     transform_size <- function(x, f_rescale) {
                       f_rescale(x) - f_rescale(0)
                     }

                     pieGrob(coords$x,
                             coords$y %||% coords$npcy,
                             r = transform_size(coords$r, panel_params$x$rescale),
                             prop = coords$prop,
                             pie_gp = grid::gpar(fill = coords$colour),
                             background_gp = grid::gpar(fill = coords$fill)
                     )

                   }

)

#'
#' @export
geom_pie <- function(mapping = NULL, data = NULL, stat = "identity",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, ...) {

  ggplot2::layer(
    geom = GeomPie, mapping = mapping,  data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}




#' @export
#' @import grid
pieGrob <- function(x, y,
                    r, prop,
                    default.units = "npc",
                    name = NULL,
                    pie_gp = gpar(),
                    background_gp = gpar(),
                    vp = NULL) {
  # Check if input needs to be converted to units
  if (!is.unit(x)) {
    x <- unit(x, default.units)
  }
  if (!is.unit(y)) {
    y <- unit(y, default.units)
  }
  if (!is.unit(r)) {
    r <- unit(r, default.units)
  }

  stopifnot(prop >= 0 & prop <= 1)
  gTree(
    x = x,
    y = y,
    r = r,
    prop = prop,
    name = name,
    pie_gp = pie_gp,
    background_gp = background_gp,
    vp = vp,
    cl = "pie"
  )
}

#' @export
makeContent.pie <- function(x) {
  x_pos <- convertWidth(x$x, unitTo = "cm", valueOnly = TRUE)
  y_pos <- convertHeight(x$y, unitTo = "cm", valueOnly = TRUE)
  r <- convertWidth(x$r, unitTo = "cm", valueOnly = TRUE)


  circle_grob <- circleGrob(x = x_pos,
                            y = y_pos,
                            r = r,
                            default.units = "cm",
                            gp = gpar(fill = x$background_gp$fill))

  pie_data <- Map(function(x, y, r, prop, id) {
    theta <- seq(0 + pi / 2,
                 - prop * 2 * pi + pi/2,
                 length.out = prop * 100)

    data.frame(
      id = id,
      x = c(x, x + r * cos(theta)),
      y = c(y, y + r * sin(theta))
    )
  },
  x = x_pos,
  y = y_pos,
  r = r,
  prop = x$prop,
  id = seq_along(x_pos))


  pie_data <- do.call(rbind, pie_data)

  pie <- polygonGrob(x = pie_data$x,
                     y = pie_data$y,
                     id = pie_data$id,
                     gp = gpar(fill = x$pie_gp$fill),
                     default.units = "cm")

  setChildren(x, gList(circle_grob, pie))
}


#' @export
StatCountPositives <- ggproto("StatCountPositives", Stat,
                              compute_group = function(data, scales) {
                                group_by(data, across(-y)) |>
                                  summarise(
                                    prop =  sum(y > 0) / n(),
                                    .groups = "drop")
                              },

                              required_aes = c("x", "y")
)

#' @export
stat_count_positives <- function(mapping = NULL, data = NULL, geom = "pie",
                                 position = "identity", na.rm = FALSE, show.legend = NA,
                                 inherit.aes = TRUE, ...) {
  layer(
    stat = StatCountPositives, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
