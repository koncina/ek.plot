#' Set absolute ggplot panel sizes (adapted from https://stackoverflow.com/a/32583612)
#'
#' @importFrom ggplot2 ggplotGrob
#' @importFrom fs dir_create
#'
#' @param x a ggplot
#' @param filename The basename of the output file
#'
#' @export
set_panel_size <- function(p = NULL, g = ggplotGrob(p),
                           width = NULL,
                           height = NULL){
  panels <- grep("panel", g$layout$name)
  panel_index_col <- unique(g$layout$l[panels])
  panel_index_row <- unique(g$layout$t[panels])
  n_col <- length(panel_index_col)
  n_row <- length(panel_index_row)



  if (!is.null(width)) {
    stopifnot(inherits(width, "unit"))
    if (n_col > 1 & length(width) == 1) width <- rep(width,  n_col)
    g$widths[panel_index_col] <- width
  }

  if (!is.null(height)) {
    stopifnot(inherits(height, "unit"))
    if (n_row > 1 & length(height) == 1) height <- rep(height,  n_row)
    g$heights[panel_index_row] <- height
  }

  g
}


#' Save a ggplot or ComplexHeatmap as a standalone file
#' and return path to the file.
#' I tries to guess the size of the plot if possible
#' (downstream of set_panel_size() or after defining a heatmap size)
#'
#' @importFrom ggplot2 ggsave ggplotGrob
#' @importFrom fs dir_create
#' @importFrom grid convertWidth convertHeight unitType
#'
#' @param x a ggplot or ComplexHeatmap.
#' @param filename The basename of the output file
#' @param device Device to use. Defaults to PNG.
#' @param width,height,units The settings for the output file.
#' @param ... Arguments passed to device or draw.
#'
#' @export
write_plot <- function (x, filename, device = grDevices::png,
                        width = NA, height = NA, units = "cm", ...) {
  UseMethod("write_plot", x)
}


#' @export
write_plot.gtable <- function(x, filename, device = grDevices::png,
                              width = NA, height = NA, units = "cm", ...) {

  if (isTRUE(is.na(width)) && !"null" %in% unitType(x$widths))
    width <- convertWidth(sum(x$widths),
                          unitTo = units, valueOnly = TRUE)

  if (isTRUE(is.na(height)) && !"null" %in% unitType(x$heights))
    height <- convertHeight(sum(x$heights),
                            unitTo = units, valueOnly = TRUE)

  ggsave(filename, plot = x, device = device, width = width, height = height, units = units, ...)
  filename
}

#' @export
write_plot.ggplot <- function(x, filename, device = grDevices::png,
                              width = NA, height = NA, units = "cm", ...) {


  x <- ggplotGrob(x)
  write_plot(x, filename, device,
             width, height, units, ...)

}

get_heatmap_size <- function(ht, units = "px", ...) {
  ht <- draw(ht, ...)
  w  <-  ComplexHeatmap:::width(ht)
  w <-  convertX(w, units, valueOnly = TRUE)
  h <-  ComplexHeatmap:::height(ht)
  h <-  convertY(h, units, valueOnly = TRUE)
  list(width = w,
       height = h)
}

#' @export
#' @param res The resolution for the output file.
write_plot.Heatmap <- function(x, filename, device = grDevices::png,
                               width = NA, height = NA, units = "cm", res = 300, ...) {

  if (is.na(width) | is.na(height)) {
    ht_size <- get_heatmap_size(x, units = units, ...)
    width <- ht_size$width
    height <- ht_size$height
  }

  device(filename, width = width, height = height, units = units, res = res)
  draw(x, ...)
  invisible(dev.off())
  filename
}
