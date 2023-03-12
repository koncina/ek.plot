#' Set absolute ggplot panel sizes (adapted from https://stackoverflow.com/a/32583612)
#'
#' @importFrom ggplot2 ggplotGrob
#' @importFrom fs dir_create
#'
#' @param x a ggplot
#' @param filename The basename of the output file
#'
#' @export
#
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


#' Save a ggplot as a standalone file using ggsave
#'
#' @importFrom ggplot2 ggsave ggplotGrob
#' @importFrom fs dir_create
#' @importFrom grid convertWidth convertHeight unitType
#'
#' @param x a ggplot
#' @param filename The basename of the output file
#'
#' @export
write_plot <- function(x, filename, device = grDevices::png, width = NA, height = NA, units = "cm", ...) {

  f <- get_report_path(filename)

  if ("ggplot" %in% class(x))
    x <- ggplotGrob(x)

  if (isTRUE(is.na(width)) && !"null" %in% unitType(x$widths))
    width <- grid::convertWidth(sum(x$widths),
                                unitTo = units, valueOnly = TRUE)

  if (isTRUE(is.na(height)) && !"null" %in% unitType(x$heights))
    height <- grid::convertHeight(sum(x$heights),
                                  unitTo = units, valueOnly = TRUE)

  ggsave(f, plot = x, device = device, width = width, height = height, units = units, ...)
  f
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

#' Save a ComplexHeatmap as a standalone png file
#'
#' @importFrom ragg agg_png
#' @importFrom fs dir_create
#'
#' @param x a ComplexHeatmap object
#' @param filename The basename of the output file
#' @param width,height,units,res The settings for the output file
#'
#' @export
write_heatmap_png <- function(x, filename, ..., width = NA, height = NA, units = "mm", res = 150) {
  stopifnot(str_detect(filename, "\\.png$"))
  f <- get_report_path(filename)

  if (is.na(width) | is.na(height)) {
    ht_size <- get_heatmap_size(x, units = units, ...)
    width <- ht_size$width
    height <- ht_size$height
  }

  agg_png(f, width = width, height = height, units = units, res = res)
  draw(x, ...)
  invisible(dev.off())
  f
}

#' Save a ComplexHeatmap as a standalone pdf file
#'
#' @importFrom grDevices pdf
#' @importFrom fs dir_create
#'
#' @param x a ComplexHeatmap object
#' @param filename The basename of the output file
#' @param width,height The settings for the output file
#'
#' @export
write_heatmap <- function(x, filename, ..., width = NA, height = NA) {
  stopifnot(str_detect(filename, "\\.pdf$"))
  f <- get_report_path(filename)

  if (is.na(width) | is.na(height)) {
    ht_size <- get_heatmap_size(x, units = "inches", ...)
    width <- ht_size$width
    height <- ht_size$height
  }

  pdf(f, width = width, height = height)
  draw(x, ...)
  dev.off()
  f
}
