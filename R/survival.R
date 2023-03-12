NULL

round_nsmall <- function(x, digits = 2) {
  format(round(x, digits), nsmall = digits)
}

#' Tidy output for `survival::coxph()`
#'
#' Wrapper around `survival::coxph()` to provide a tidy tibble using broom.
#'
#' @importFrom survival Surv coxph cox.zph
#' @importFrom broom tidy
#'
#' @param data A tibble
#' @param time,event,predictor the variable names to peroform the COX proportional hazard model.
#' @param exponentiate Boolean telling whether to exponentiate the hazard ratio.
#' @param cox_zph Perform the `surival::cox.zph` test to check the assumption for the Cox regression model fit.
#'
#' @export
tidy_coxph <- function (data, time, event, predictor, exponentiate = TRUE, cox_zph = TRUE, ...) {
  UseMethod("tidy_coxph", data)
}


#' @export
tidy_coxph.data.frame <- function(data, time, event, predictor, exponentiate = TRUE, cox_zph = TRUE) {
  stopifnot(!is_grouped_df(data))

  surv <-Surv(pull(data, {{time}}),
              pull(data, {{event}}))
  x <- rename(data, .predictor = {{predictor}})

  n_stratus <- length(unique(pull(x, .predictor)))
  if (n_stratus > 2) {
    stop("The function does not handle more than 2 predictor levels yet", call. = FALSE)
  }

  n_by_stratus <- count(x, .predictor)
  n_by_stratus <- pivot_wider(n_by_stratus, names_from = .predictor,
                              values_from = n, names_prefix = "n_")
  coxph_fit <- with(x, coxph(surv ~ .predictor))

  x <- tidy(coxph_fit, exponentiate = {{exponentiate}}, conf.int = TRUE)
  x <- mutate(x,
              coxph = list(coxph_fit),
              across(term, str_replace,
                     ".predictor", paste0(quo_name(enquo(predictor)), "_")))
  x <- bind_cols(x, n_by_stratus)

  if (isTRUE(cox_zph)) {
    x <- mutate(x, p_val_coxzph = cox.zph(coxph_fit)[["table"]][[".predictor", "p"]])
  }

  select(x, coxph, starts_with("n_"), hr = estimate, hr_low = conf.low, hr_high = conf.high, p_val = p.value, p_val_coxzph)
}

#' @export
tidy_coxph.grouped_df <- function(data, time, event, predictor, exponentiate = TRUE, cox_zph = TRUE) {
  bind_cols(group_keys(data), map_dfr(group_split(data), tidy_coxph, {{time}}, {{event}}, {{predictor}}, {{exponentiate}}))
}
