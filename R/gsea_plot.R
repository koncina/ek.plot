plot_gsea_es <- function(data, title, palette) {

  ggplot(data,
         aes(x = x,
             y = running_score,
             colour = description)) +
    geom_hline(yintercept = 0, colour = "darkgray") +
    geom_line() +
    theme_linedraw(10) +
    theme(axis.text.y = element_markdown(),
          legend.position = "none",
          axis.title.x = element_blank(),
          plot.margin = margin(5, 5, 0, 5),
          plot.title = element_markdown(size = 10, hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(colour = "gray", size = 0.2),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), limits = c(0, max(data$x))) +
    scale_y_continuous(limits = \(x) c(-max(abs(x)), max(abs(x))),
                       minor_breaks = c(-1, -0.5, 0.5, 1),
                       breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                       labels = c("-1", "", "-0.5", "", "0", "", "0.5", "", "1")) +
    scale_colour_manual(values = palette) +
    labs(y = "Enrichment score",
         title = {{title}},
         colour = NULL)
}

plot_gsea_ranks <- function(data, palette) {

  n_genes <- max(data$x)

  data <-  filter(data, position != 0) |>
    arrange(gs_id) |>
    mutate(
      lbl = gs_lbl,#glue::glue("{gs_lbl} ({gs_id})"),
      across(gs_id, as.numeric),
      gs_id = 1 + max(gs_id) - gs_id,
      across(lbl, compose(fct_rev, fct_inorder))
    )

  ggplot(data,
         aes(x = x,
             y = lbl,
             colour = description)) +
    geom_blank() +
    geom_segment(aes(xend = x, y = gs_id - 0.5, yend = gs_id + 0.5),
                 linewidth = 0.2) +
    scale_y_discrete(expand = expansion(add = 0.5), position = "right") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                       limits = c(0, n_genes),
                       labels = \(x) x / 1000) +
    scale_colour_manual(values = palette) +
    theme_linedraw(10) +
    theme(axis.ticks.y = element_blank(),
          legend.position = "none",

          axis.title.y = element_blank(),
          axis.text.y.right = element_markdown(),
          plot.margin = margin(0, 5, 5, 5),
          axis.title.x = element_markdown(),
          panel.grid.major.y = element_blank())  +
    labs(x = "Gene rank (x10<sup>3</sup>)")
}

plot_gsea_gradient <- function(data, lbl) {

  distinct(data, x) |>
    ggplot(aes(x = x, y = 1, fill = -x, colour = -x)) +
    geom_tile(stroke = 0) +
    scale_colour_gradientn(colours = c("#377EB8", "white", "#E41A1C")) +
    theme_linedraw(10) +
    theme(axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(0, 5, 0, 5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    scale_y_continuous(expand = expansion()) +
    annotate("richtext", x = 1, y = 1, label = lbl[1], hjust = 0,
             fill = NA, label.color = NA,
             label.margin = unit(c(0, 0, -1, 1), "mm"),
             family = "sans",
             size = 8 / ggplot2:::.pt) +
    annotate("richtext", x = max(data$x), y = 1, label = lbl[2], hjust = 1,
             fill = NA, label.color = NA,
             label.margin = unit(c(0, 1, -1, 0), "mm"),
             family = "sans",
             size = 8 / ggplot2:::.pt) +
    scale_x_continuous(expand = expansion(), limits = c(0, max(data$x))) +
    labs(x = "Gene rank")
}


plot_gsea_volcano <- function(data, palette) {

  max_x <- pull(data, nes) |>
    abs() |>
    max()

  ggplot(data, aes(x = nes, y = -log10(p_adjust))) +
    geom_point(aes(colour = gs_id, alpha = !is.na(gs_id)), size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits = c(-max_x, max_x), breaks = -3:3) + #
    scale_y_continuous(
      limits = c(0, 10),
      expand = expansion(add = c(0, 1)),
      oob = scales::squish,
      breaks = 0:5*2
    ) +
    scale_colour_manual(values = palette) + #c(`FALSE` = "black", `TRUE` = "red")) +
    scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 1)) +
    theme_classic(10) +
    theme(legend.position = "none",
          axis.title.y = element_markdown(margin = margin()),
          plot.title = element_markdown(),
          plot.background = element_rect(fill = NA, colour = NA),
          #plot.margin = margin(t = 15),
          plot.margin = margin(l = unit(2, "mm"), r = unit(5, "pt"))

          ) +
    labs(x = "NES",
         y = "-log<sub>10</sub>(p-adj)")
}

#' Prepare the GSEA plot
#'
#' @importFrom tibble deframe
#' @importFrom dplyr filter pull
#' @importFrom purrr map_dfr
#' @importFrom enrichplot gsInfo
#' @importFrom janitor make_clean_names
#' @import ggplot2
#'
#' @param gsea,gs_tidy The GSEA and selected genesets to be highlighted
#' @param title Title of the plot
#' @param lbl vector of size 2 with labels to be shown on the gradient bar (DGE contrast)
#'
#' @export
plot_gsea <- function(
    gsea, gs_tidy,
    title = "Gene set",
    lbl = c("Left", "Right"),
    palette) {

  palette <- select(gs_tidy, gs_name, colour) |> deframe()

  data_summary <- as_tibble(gsea,
                            .name_repair = janitor::make_clean_names) |>
    left_join(gs_tidy, by = c("description" = "gs_name")) |>
    arrange(desc(nes), p_adjust) |>
    mutate(#gs_id = replace(NA_integer_, !is.na(gs_lbl),
      #               seq_along(nes[!is.na(gs_lbl)])),
      gs_id =  replace(NA_integer_, !is.na(gs_lbl), description[!is.na(gs_lbl)]),
      #across(gs_id, as.factor)
      across(gs_id, fct_inorder)
    )


  data <- filter(gs_tidy, !is.na(gs_lbl)) |>
    pull(gs_name) |>
    map_dfr(\(x) enrichplot:::gsInfo(object = gsea, geneSetID = x)) |>
    as_tibble(.name_repair =  janitor::make_clean_names) |>
    left_join(data_summary,
              by = "description") |>
    arrange(gs_id)

  list(
    volcano = plot_gsea_volcano(data_summary, palette),
    es = plot_gsea_es(data, title, palette) /
      plot_gsea_gradient(data, lbl = lbl) /
      plot_gsea_ranks(data, palette),
    genesets = gs_tidy,
    data_summary,
    data
  )
}

#' Build the GSEA montage using the output of `plot_gsea_volcano()`:
#' Showing the volcano plot on upper right corner
#' And a montage for selected genesets showing the running enrichment score,
#' gradient and ranks.
#'
#' @importFrom ggplot2 ggplotGrob
#' @importFrom patchwork patchworkGrob
#' @importFrom dplyr n_distinct
#' @importFrom grid convertWidth convertHeight unitType rectGrob
#' @importFrom purrr map_dbl
#' @importFrom gtable gtable gtable_add_grob
#'
#' @param p a GSEA plot generated by `plot_gsea_volcano()`
#'
#' @export
build_gsea_montage <- function(p) {
  g1 <- patchworkGrob(p[[2]]) |>
    set_panel_size(
      p = NULL,
      width = unit(3, "cm"),
      height = unit(
        c(2.6, 0.75,
          n_distinct(p[[3]]$gs_lbl) * 0.3),
        c("cm",
          "lines",
          rep("cm", n_distinct(p[[3]]$gs_lbl)))
      )
    )

  col <- unique(g1$layout$l[panels])
  row <- unique(g1$layout$t[panels])
  stopifnot(length(col) == 1 & col > 1)

  g1_x <- list(1:(col - 1), col, -(1:col)) |>
    map_dbl(
      \(x) convertWidth(
        g1$widths[x],
        unitTo = "cm",
        valueOnly = TRUE
      ) |>
        sum()
    )

  # Volcano plot (g2)

  g2 <- p[["volcano"]] |>
    set_panel_size(width = unit(2, "cm"))

  g2_w <- convertWidth(
    g2$widths,
    unitTo = "cm",
    valueOnly = TRUE
  ) |>
    sum()

  if (g2_w > g1_x[3]) {
    g1_x <- c(g1_x, g2_w - g1_x[3])
    g1_r <- 3
    g2_r <- 4
  } else {
    g1_x <- c(g1_x[1:2], g2_w, g1_x[3] - g2_w)
    g1_r <- 4
    g2_r <- 3
  }

  g1_y <- list(1:(row[1] - 1), row[1], (row[1] + 1):row[2], -(1:row[2])) |>
    map_dbl(
      \(x) convertHeight(
        g1$heights[x],
        unitTo = "cm",
        valueOnly = TRUE
      ) |>
        sum()
    )

  g2 <- ggplotGrob(p[[1]])

  g2_h <- convertHeight(
    g2$heights,
    unitTo = "cm",
    valueOnly = TRUE
  ) |> sum()

  g2 <- set_panel_size(g = g2, width = unit(2, "cm"), height = unit(3 - g2_h, "cm"))

  gt <- gtable(
    widths = unit(g1_x, "cm"),
    heights = unit(g1_y, "cm")
  )

  background <- rectGrob(gp = gpar(fill = "white", col = "white"))

  if (g1_r == 3) {
    gt <- gtable_add_grob(gt, background, t = 1, b = 4, l = 1, r = 4)
  }

  gt <- gtable_add_grob(gt, g1, t = 1, b = 4, l = 1, r = g1_r)
  gt <- gtable_add_grob(gt, g2, t = 2, b = 3, l = 3, r = g2_r)

  gt
}

