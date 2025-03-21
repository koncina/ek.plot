NULL

#' @importFrom ggplot2 aes ggproto Stat
#' @import ggplot2
#' @importFrom scales pvalue_format
#' @rdname stat_cox
#' @export
StatCox <- ggplot2::ggproto("StatCox", ggplot2::Stat,

                            compute_panel = function(data, scales) {
                              data
                              mutate(data, label = NA)
                            },

                            finish_layer = function(data, params) {
                              colour_key <- distinct(data, group, colour) |> deframe()
                              out <- group_by(data, PANEL) |>
                                tidy_coxph(time, status, group) |>
                                group_by(PANEL) |>
                                # Temporary workaround in case groups differ in each panel
                                # This will be fixed upstream (tidy_coxph())
                                pivot_longer(names_to = "g", values_to = "n", starts_with("n_"), names_prefix = "n_", values_drop_na = TRUE) |>
                                mutate(group_idx = seq_along(g)) |>
                                pivot_wider(names_from = "group_idx", values_from = c("n", "g")) |>
                                #group_by(PANEL) |>
                                mutate(
                                  p_val = scales::pvalue_format(add_p = TRUE)(p_val),
                                  hr = paste0("HR=", round(hr, 2)),
                                  n = glue::glue(
                                    paste0(
                                      "n=<span style = 'color:",
                                      colour_key[g_1],";'>{", n_1, "}</span> / <span style = 'color:",
                                      colour_key[g_2],";'>{", n_2, "}</span>"
                                    )
                                  ),
                                  # Hardcoded for now...
                                  label = paste0(hr, ", ", p_val, "<BR>", n)
                                  #label = paste(hr, p_val, n, sep = "<BR>")#,
                                  #x = 100, y = 0.5
                                ) |>
                                select(PANEL, label)

                              out <- select(data, -time, -status, -group, -colour, -label) |>
                                distinct() |>
                                left_join(out, by = "PANEL") |>
                                mutate(group = 1, colour = "black")
                              out

                            },


                            default_aes = ggplot2::aes(label = after_stat(label), x = 0, y = 0),

                            required_aes = c("time", "status")
)


stat_cox <- function(mapping = NULL, data = NULL, geom = "richtext",
                     show.legend = FALSE, inherit.aes = TRUE) {
  ggplot2::layer(
    stat = StatCox,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = FALSE,
    inherit.aes = inherit.aes,
    params = list(...)
  )

}
