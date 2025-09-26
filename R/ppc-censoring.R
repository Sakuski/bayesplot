#' PPC censoring
#'
#' @description Compare the empirical distribution of censored data `y` to the
#'   distributions of simulated/replicated data `yrep` from the posterior
#'   predictive distribution. See the **Plot Descriptions** section, below, for
#'   details.
#'
#'   Although some of the other \pkg{bayesplot} plots can be used with censored
#'   data, `ppc_km_overlay()` is currently the only plotting function designed
#'   *specifically* for censored data. We encourage you to suggest or contribute
#'   additional plots at
#'   [github.com/stan-dev/bayesplot](https://github.com/stan-dev/bayesplot).
#'
#' @name PPC-censoring
#' @family PPCs
#'
#' @template args-y-yrep
#' @param size,alpha Passed to the appropriate geom to control the appearance of
#'   the `yrep` distributions.
#' @param ... Currently only used internally.
#'
#' @template return-ggplot
#'
#' @section Plot Descriptions:
#' \describe{
#'   \item{`ppc_km_overlay()`}{
#'   Empirical CCDF estimates of each dataset (row) in `yrep` are overlaid, with
#'   the Kaplan-Meier estimate (Kaplan and Meier, 1958) for `y` itself on top
#'   (and in a darker shade). This is a PPC suitable for right-censored `y`.
#'   Note that the replicated data from `yrep` is assumed to be uncensored. Left
#'   truncation (delayed entry) times for `y` can be specified using
#'   `left_truncation_y`.
#'   }
#'   \item{`ppc_km_overlay_grouped()`}{
#'    The same as `ppc_km_overlay()`, but with separate facets by `group`.
#'   }
#' }
#'
#' @templateVar bdaRef (Ch. 6)
#' @template reference-bda
#' @template reference-km
#'
#' @examples
#' \donttest{
#' color_scheme_set("brightblue")
#'
#' # For illustrative purposes, (right-)censor values y > 110:
#' y <- example_y_data()
#' status_y <- as.numeric(y <= 110)
#' y <- pmin(y, 110)
#'
#' # In reality, the replicated data (yrep) would be obtained from a
#' # model which takes the censoring of y properly into account. Here,
#' # for illustrative purposes, we simply use example_yrep_draws():
#' yrep <- example_yrep_draws()
#' dim(yrep)
#'
#' # Overlay 25 curves
#' ppc_km_overlay(y, yrep[1:25, ], status_y = status_y)
#'
#' # With extrapolation_factor = 1 (no extrapolation)
#' ppc_km_overlay(y, yrep[1:25, ], status_y = status_y, extrapolation_factor = 1)
#'
#' # With extrapolation_factor = Inf (show all posterior predictive draws)
#' ppc_km_overlay(y, yrep[1:25, ], status_y = status_y, extrapolation_factor = Inf)
#'
#' # With separate facets by group:
#' group <- example_group_data()
#' ppc_km_overlay_grouped(y, yrep[1:25, ], group = group, status_y = status_y)
#'
#' # With left-truncation (delayed entry) times:
#' min_vals <- pmin(y, apply(yrep, 2, min))
#' left_truncation_y <- rep(0, length(y))
#' condition <- y > mean(y) / 2
#' left_truncation_y[condition] <- pmin(
#'   runif(sum(condition), min = 0.6, max = 0.99) * y[condition],
#'   min_vals[condition] - 0.001
#' )
#' ppc_km_overlay(y, yrep[1:25, ], status_y = status_y,
#'               left_truncation_y = left_truncation_y)
#' }
NULL

#' @export
#' @rdname PPC-censoring
#' @param status_y The status indicator for the observations from `y`. This must
#'   be a numeric vector of the same length as `y` with values in \{0, 1\} (0 =
#'   right censored, 1 = event).
#' @param left_truncation_y Optional parameter that specifies left-truncation
#'   (delayed entry) times for the observations from `y`. This must be a numeric
#'   vector of the same length as `y`. If `NULL` (default), no left-truncation
#'   is assumed.
#' @param extrapolation_factor A numeric value (>=1) that controls how far the
#'   plot is extended beyond the largest observed value in `y`. The default
#'   value is 1.2, which corresponds to 20 % extrapolation. Note that all
#'   posterior predictive draws may not be shown by default because of the
#'   controlled extrapolation. To display all posterior predictive draws, set
#'   `extrapolation_factor = Inf`.
#'
ppc_km_overlay <- function(
    y,
    yrep,
    yimp = NULL,
    ...,
    status_y,
    left_truncation_y = NULL,
    extrapolation_factor = 1.2,
    size = 0.25,
    alpha = 0.7
) {
  # 1. Argument and package checks (mostly from original function)
  # =============================================================
  check_ignored_arguments(..., ok_args = "add_group")
  add_group <- list(...)$add_group

  suggested_package("survival")
  suggested_package("ggfortify")
  suggested_package("dplyr")
  suggested_package("tidyr")
  suggested_package("tibble")


  if (!is.numeric(status_y) || length(status_y) != length(y) || !all(status_y %in% c(0, 1))) {
    stop("`status_y` must be a numeric vector of 0s and 1s the same length as `y`.", call. = FALSE)
  }
  if (!is.null(yimp) && (!is.numeric(yimp) || length(yimp) != length(y))) {
    stop("`yimp` must be a numeric vector the same length as `y`.", call. = FALSE)
  }
  if (!is.null(left_truncation_y)) {
    if (!is.numeric(left_truncation_y) || length(left_truncation_y) != length(y)) {
      stop("`left_truncation_y` must be a numeric vector of the same length as `y`.", call. = FALSE)
    }
  }
  if (extrapolation_factor < 1) {
    stop("`extrapolation_factor` must be greater than or equal to 1.", call. = FALSE)
  }
  if (extrapolation_factor == 1.2) {
    message(
      "Note: `extrapolation_factor` now defaults to 1.2 (20%).\n",
      "To display all posterior predictive draws, set `extrapolation_factor = Inf`."
    )
  }

  # 2. Manual Data Preparation (replaces ppc_data())
  # =================================================
  # Create individual tibbles for y, yimp, and yrep
  y_df <- tibble::tibble(
    value = y,
    y_id = seq_along(y),
    source = "y",
    rep_id = 0,
    event = status_y # Use the actual status for observed data
  )

  yimp_df <- NULL
  if (!is.null(yimp)) {
    yimp_df <- tibble::tibble(
      value = yimp,
      y_id = seq_along(yimp),
      source = "yimp",
      rep_id = 0,
      event = 1 # Imputed data is treated as a complete event
    )
  }

  S <- nrow(yrep)
  N <- ncol(yrep)
  yrep_df <- tibble::as_tibble(yrep, .name_repair = ~ as.character(1:N)) %>%
    dplyr::mutate(rep_id = 1:S) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of("rep_id"),
      names_to = "y_id",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      y_id = as.integer(.data$y_id),
      source = "yrep",
      event = 1 # Replicated data is treated as a complete event
    )

  # Combine into a single long-format data frame
  combined_data <- dplyr::bind_rows(y_df, yimp_df, yrep_df) %>%
    dplyr::mutate(
      # Create a unique identifier for each curve for the survfit formula
      strata_group = dplyr::case_when(
        .data$source == "y"    ~ "italic(y)",
        .data$source == "yimp" ~ "italic(y^imp)",
        .data$source == "yrep" ~ paste0("[rep] (", .data$rep_id, ")")
      )
    )

  # 3. Kaplan-Meier Estimation
  # ============================
  if (is.null(left_truncation_y)) {
    sf_form <- survival::Surv(time = combined_data$value, event = combined_data$event) ~ strata_group
  } else {
    sf_form <- survival::Surv(time = left_truncation_y[combined_data$y_id], time2 = combined_data$value, event = combined_data$event) ~ strata_group
  }

  if (!is.null(add_group)) {
    add_group_df <- tibble::tibble(y_id = seq_along(y), add_group = add_group)
    combined_data <- dplyr::inner_join(combined_data, add_group_df, by = "y_id")
    sf_form <- update(sf_form, . ~ . + add_group)
  }

  sf <- survival::survfit(sf_form, data = combined_data)

  # Rename strata for ggfortify parsing if add_group is used
  if (!is.null(add_group)) {
    names(sf$strata) <- sub("add_group=", "add_group:", names(sf$strata))
  }

  # 4. Prepare fortified data for ggplot
  # ======================================
  fsf <- fortify(sf)

  # Handle splitting strata if add_group was used
  if(any(grepl("add_group", levels(fsf$strata)))){
    strata_split <- strsplit(as.character(fsf$strata), split = ", add_group:")
    fsf$strata <- as.factor(sapply(strata_split, "[[", 1))
    fsf$group <- as.factor(sapply(strata_split, "[[", 2))
  }

  # Create a column to map aesthetics (color, size, alpha)
  fsf <- fsf %>%
    dplyr::mutate(
      curve_type = dplyr::case_when(
        grepl("italic\\(y\\)$", .data$strata) ~ "y",
        grepl("y\\^imp", .data$strata)      ~ "yimp",
        TRUE                               ~ "yrep"
      )
    )

  # Filter yrep curves based on the extrapolation factor
  max_time_y <- max(y, na.rm = TRUE)
  fsf <- fsf %>%
    dplyr::filter(.data$curve_type == "y" | .data$time <= max_time_y * extrapolation_factor)

  # Ensure correct plotting order: y on top, then yimp, then yrep at the bottom
  y_level <- "italic(y)"
  yimp_level <- if (!is.null(yimp)) "italic(y^imp)" else NULL
  yrep_levels <- unique(fsf$strata[fsf$curve_type == "yrep"])
  new_levels <- c(yrep_levels, yimp_level, y_level)
  fsf$strata <- factor(fsf$strata, levels = new_levels)

  # 5. Generate the Plot using multiple layers for correct ordering
  # ===============================================================
  p <- ggplot(
    data = fsf,
    mapping = aes(
      x = .data$time,
      y = .data$surv,
      group = .data$strata,
      color = .data$curve_type,
      size = .data$curve_type,
      alpha = .data$curve_type
    )
  ) +
    # Layer 1 (Bottom): yrep curves, filtered using the dot-pipe pronoun
    geom_step(data = . %>% dplyr::filter(curve_type == "yrep"), na.rm = TRUE) +
    # Layer 2 (Middle): yimp curve
    geom_step(data = . %>% dplyr::filter(curve_type == "yimp"), na.rm = TRUE) +
    # Layer 3 (Top): y curve
    geom_step(data = . %>% dplyr::filter(curve_type == "y"), na.rm = TRUE) +
    # Reference lines
    hline_at(0.5, linewidth = 0.1, linetype = 2, color = get_color("dh")) +
    hline_at(c(0, 1), linewidth = 0.2, linetype = 2, color = get_color("dh")) +
    # Manual scales to control the appearance of each curve type
    scale_color_manual(
      name = NULL,
      values = c("y" = get_color("dark"), "yimp" = get_color("mid"), "yrep" = get_color("light")),
      labels = c("y" = expression(italic(y)), "yimp" = expression(italic(y)^{imp}), "yrep" = expression(italic(y)^{rep}))
    ) +
    scale_size_manual(
      name = NULL,
      values = c("y" = 1, "yimp" = 1, "yrep" = size),
      guide = "none"
    ) +
    scale_alpha_manual(
      name = NULL,
      values = c("y" = 1, "yimp" = 0.75, "yrep" = alpha),
      guide = "none"
    ) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    # Standard bayesplot theme elements
    xlab(y_label()) +
    yaxis_title(FALSE) +
    xaxis_title(FALSE) +
    yaxis_ticks(FALSE) +
    bayesplot_theme_get() +
    theme(legend.text.align = 0)

  return(p)
}

#' @export
#' @rdname PPC-censoring
#' @template args-group
ppc_km_overlay_grouped <- function(
  y,
  yrep,
  group,
  ...,
  status_y,
  left_truncation_y = NULL,
  extrapolation_factor = 1.2,
  size = 0.25,
  alpha = 0.7
) {
  check_ignored_arguments(...)

  p_overlay <- ppc_km_overlay(
    y = y,
    yrep = yrep,
    add_group = group,
    ...,
    status_y = status_y,
    left_truncation_y = left_truncation_y,
    size = size,
    alpha = alpha,
    extrapolation_factor = extrapolation_factor
  )

  p_overlay +
    facet_wrap("group") +
    force_axes_in_facets()
}
