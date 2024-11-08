source("helper-functions.R")

fce_ind_release <- function(data, release_sched, mark_col = pit, ...) {
  data |>
    filter(
      !is.na({{ mark_col }}),
      ...
    ) |>
    left_join(release_sched, by = join_by(date == collect_date)) |>
    select(release_date, species, lifestage, {{ mark_col }}, count)
}

fce_ind_release_summ <- function(release, trap_op) {
  release |>
    group_by(release_date, species, lifestage) |>
    summarize(
      count = sum(count),
      .groups = "drop"
    ) |>
    left_join(trap_op, by = join_by(release_date == date)) |>
    arrange(release_date, species, lifestage) |>
    select(release_date, species, lifestage, count, doy) |>
    mutate(release_idx = seq_along(release_date))
}

fce_ind_recap <- function(data, mark_col = pit, ...) {
  data |>
    filter(
      !is.na({{ mark_col }}),
      ...
    ) |>
    select(recapture_date = date, {{ mark_col }})
}

fce_unmarked <- function(data, trap_op, ...) {
  cdat <- data |>
    filter(...) |>
    select(date, species, lifestage, count) |>
    summarize(
      count = sum(count),
      .by = c(date, species, lifestage)
    )
  left_join(trap_op, cdat, by = join_by(date))
}

fce_release_df <- function(release_ind) {
  release_ind |>
    summarize(
      count = sum(count),
      .by = c(release_date, species, lifestage)
    ) |>
    ## left_join(recap_window, by = join_by(species)) |>
    mutate( ## rec_idx_start = cumsum(recap_window + 1) - recap_window,
      ## Don't need to subtract one because we need an extra spot for the
      ## lost count
      ## rec_idx_end = rec_idx_start + recap_window,
      release_idx = seq_along(release_date)
    )
}

## Find the first and last day that the trap operates each year
fce_trap_season <- function(data) {
  data |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarize(
      start = min(date),
      end = max(date)
    )
}

## Extract the operational history of the trap.
fce_trap_op <- function(data, ...) {
  ## Extract the trap season for each year. May be different lengths per year,
  ## will get the row indices later.
  trap_season <- fce_trap_season(data)

  ## Assume that collection is occurring any day when a row is present in the
  ## data set
  trap_coll <- data |>
    summarize(op = TRUE, .by = date)

  trap_season |>
    mutate(date = map2(start, end, seq, by = "1 day")) |>
    select(year, date) |>
    unnest(date) |>
    left_join(trap_coll, by = join_by(date)) |>
    mutate(
      op = replace_na(op, FALSE),
      doy = yday(date)
    ) |>
    filter(...)
}

## Generate the data frame with individual mark-recapture histories
fce_ind <- function(release_ind, recap_ind, mark_col = pit) {
  left_join(release_ind, recap_ind, join_by({{ mark_col }})) |>
    select(release_date, recapture_date, species, lifestage, pit, count) |>
    mutate(travel_time = recapture_date - release_date)
}

## Define the species-specific default recapture windows so that it doesn't have
## to be passed around a bunch.
fce_recap_window_default <- function(release_ind, recap_ind) {
  fce_release_df(release_ind) |>
    select(release_date, species) |>
    left_join(
      tribble(
        ~species, ~recap_window,
        "Steelhead", 35,
        "Coho", 35,
        "Chinook", 105
      ),
      by = join_by(species)
    )
}

## Define the recapture window based on the maximum observed recapture delay
fce_recap_window_obs <- function(release_ind, recap_ind) {
  rw_df <- left_join(release_ind, recap_ind, by = join_by(pit)) |>
    filter(!is.na(recapture_date)) |>
    mutate(recap_window = recapture_date - release_date) |>
    summarize(
      recap_window = max(recap_window),
      .by = c(species, lifestage)
    ) |>
    mutate(recap_window = as.numeric(recap_window) + 1)

  fce_release_df(release_ind) |>
    select(release_date, species, lifestage) |>
    left_join(rw_df, by = join_by(species, lifestage))
}

## Changes multinomial observation likelihood too much; shouldn't use
## fce_recap_window_ragged <- function(release_ind, recap_ind) {
##   fce_ind(release_ind, recap_ind) |>
##     filter(!is.na(recapture_date)) |>
##     summarize(recap_window = as.numeric(max(travel_time)),
##               .by = c(release_date, species))

## }

fce_ragged_index <- function(lengths, additional = 0) {
  cumsum(c(1, lengths + additional))
}

fce_recap_df <- function(release_ind, recap_ind,
                         trap_op,
                         recap_window_fun = fce_recap_window_obs,
                         drop_notrap = TRUE) {
  ## Summarize release information. Probably already did this outside the
  ## function, but need the individual histories so it's not a big deal to do it
  ## again. Also attach the recapture window for each release.
  release_df <- fce_release_df(release_ind) |>
    left_join(recap_window_fun(release_ind, recap_ind),
      by = join_by(release_date, species, lifestage)
    )

  ## Create the skeleton for recapture information; fill in the potential dates
  ## of recapture based on the release date and the release window
  recap_df <- release_df |>
    mutate(recapture_date = map2(
      release_date, recap_window,
      \(rd, rw) c(rd + 1:rw)
    )) |>
    select(release_date, species, lifestage, recapture_date) |>
    unnest(recapture_date) |>
    ## Filter out recaptures that occur after the end of the trapping season
    filter(recapture_date %in% trap_op$date) |>
    mutate(
      doy = yday(recapture_date),
      year = year(recapture_date),
      travel_time = as.numeric(recapture_date - release_date)
    )

  ## Summarize the recapture information as number recaptured each recapture day
  recap_summ <- left_join(release_ind, recap_ind, by = join_by(pit)) |>
    summarize(
      count = sum(count),
      .by = c(species, lifestage, release_date, recapture_date)
    )

  ## Associate the recapture counts by day with the specific days, filling in
  ## zeros where necessary
  rdf <- left_join(recap_df, recap_summ,
    by = join_by(
      species,
      lifestage,
      release_date,
      recapture_date
    )
  ) |>
    mutate(count = replace_na(count, 0)) |>
    left_join(select(trap_op, date, op),
      by = join_by(recapture_date == date)
    )

  ## Add option to drop rows (and thus any likelihood calculations) of rows
  ## where the trap is closed and thus there is probability zero of any
  ## captures. The zeros add a constant to the posterior so they don't need to
  ## be included. This also allows probabilities to be calculated in logit space
  ## which could provide an improvement in numerical stability. This *is*
  ## accounted for when calculating the probabiltiy of recapture in the model
  ## (where the majority of the p_recap vector is calculated, currently line
  ## 183).
  if (drop_notrap) {
    rdf <- rdf |>
      filter(op)
  }
  return(rdf)
}

fce_num_lost <- function(release_df, recap_df) {
  n_recap <- recap_df |>
    summarize(
      recap_count = sum(count),
      .by = c(release_date, species, lifestage)
    )
  left_join(release_df, n_recap,
    by = join_by(release_date, species, lifestage)
  ) |>
    mutate(
      n_lost = count - recap_count
    ) |>
    pluck("n_lost")
}

## fce_recap_index <- function(recap_df) {
##   recap_df |>
##     mutate(index = seq_along(recapture_date)) |>
##     group_by(release_date, species) |>
##     summarize(rec_idx_start = min(index),
##               rec_idx_end = max(index),
##               .groups = "drop")
## }

## fce_recap_window_df <- function(trap_op, rel_day, recap_window, list = TRUE) {
##   rw <- trap_op |>
##     arrange(date) |>
##     filter(date > rel_day) |>
##     slice_head(n = recap_window) |>
##     rename(recapture_date = date) |>
##     bind_rows(tibble(recapture_date = NA, doy = NA, op = NA)) |>
##     mutate(travel_time = seq_along(recapture_date))
##   if (list) {
##     rw <- list(rw)
##   }
##   rw
## }

fce_knot_vec <- function(inner, outer, expand = 1.1) {
  inr <- range(inner)
  inr <- inr + c(-1, 1) * diff(inr) * (expand - 1)
  outr <- range(outer)
  outr <- outr + c(-1, 1) * diff(outr) * (expand - 1)
  c(outr[1], inr, outr[2])
}

fce_parse_formula <- function(formula, data, knots = NULL) {
  if (attr(terms(formula), "response")) {
    warning("Must use a one-sided formula, response removed")
    formula <- formula[-2]
  }
  gp <- interpret.gam(formula)
  sm <- lapply(gp$smooth.spec,
    smoothCon,
    data = data,
    knots = knots,
    absorb.cons = TRUE,
    scale.penalty = FALSE,
    diagonal.penalty = TRUE
  )
  list(
    linear = gp$pf,
    smooth = sm,
    data = data,
    pred_names = gp$pred.names
  )
}

fce_construct_mm <- function(pform, data = NULL,
                             check_fullrank = FALSE) {
  if (is.null(data)) {
    data <- pform$data
  } else {
    ## Make sure that all required columns are present
    if (!all(pform$pred_names %in% names(data))) {
      mc <- setdiff(pform$pred_names, names(data))
      stop("data is missing column(s): ", mc)
    }
    ## Makes sure that any columns that are factors in the original data set are
    ## also factors with the same set of levels in the new data.
    fct_cols <- names(Filter(is.factor, pform$data))
    data <- data |>
      mutate(
        across(
          all_of(fct_cols),
          \(col) factor(col, levels = levels(pform$data[[cur_column()]]))
        )
      )
  }

  lin_mm <- model.matrix(pform$linear, data = data, na.action = na.fail)
  ## The intercept in attr(terms(pform$linear)) has a numeric value of zero, so
  ## concatenate it at the start of the vector and then add one to the indices
  ## used below
  lin_labs <- c("Intercept", attr(terms(pform$linear), "term.labels")) |>
    set_names()
  lin_scl <- lin_labs[attr(lin_mm, "assign") + 1]

  sm <- flatten(pform$smooth)
  sm_lab <- map_chr(sm, pluck, "label")
  ## sm <- set_names(sm, sm_lab)
  sm_term <- map(sm, pluck, "term")
  sm_term <- map_chr(sm_term, paste0, collapse = ":")
  sm_pnull <- map(sm, \(s) c(rep(FALSE, s$rank), rep(TRUE, s$null.space.dim)))

  sm_scl <- map2(
    sm_term, sm_pnull,
    \(t, p) paste0(t, ifelse(p, "_null", ""))
  )

  sm_mm <- map(sm, PredictMat, data = data)
  sm_mm <- map2(
    sm_mm, sm_lab,
    \(s, l) {
      colnames(s) <- rep(l, ncol(s))
      s
    }
  )

  ## Need `.init = vector()` for the case where there are no smoothers to
  ## `cbind`
  mm <- cbind(lin_mm, reduce(sm_mm, cbind, .init = vector()))

  if (check_fullrank && Matrix::rankMatrix(mm) < ncol(mm)) {
    warning("Model matrix is not full rank")
  }

  scl_idx <- flatten_chr(c(lin_scl, sm_scl))
  ## scl_idx <- factor(scl_idx, levels = unique(scl_idx))
  attr(mm, "scale_index") <- scl_idx

  mm
}

## Create the mgcv smoother objects that we can then use to get design matrices
fce_spl <- function(spec, data, knots = NULL) {
  smoothCon(
    spec,
    data = data,
    knots = knots,
    diagonal.penalty = TRUE, # Don't need to pass penalty matrix to Stan
    scale.penalty = FALSE
  )
}

## Get the matrices
fce_basis <- function(spl, data) {
  basis <- map(spl, PredictMat, data = data)
  Reduce(cbind, basis)
}

fce_runsize_basis <- function(spec, data, knots = NULL) {
  spl <- smoothCon(
    spec,
    data = data,
    diagonal.penalty = TRUE,
    scale.penalty = FALSE
  )
  basis <- map(spl, PredictMat, data = data)
  Reduce(cbind, basis)
}

get_travel_form_idx <- function(travel_form) {
  form <- switch(travel_form,
    "log-normal" = ,
    "lognormal" = ,
    "log_normal" = ,
    "log normal" = 1L,
    "inverse-gaussian" = ,
    "inv-gaussian" = ,
    "inv-gauss" = ,
    "inverse_gaussian" = ,
    "inv_gaussian" = ,
    "inv_gauss" = ,
    "inverse gaussian" = 2L,
    "exponential" = ,
    "exp" = 3L
  )
  npar <- switch(form,
    2, ## log-normal,
    2, ## inverse-gaussian
    1
  ) ## exponential
  list(form = form, npar = npar)
}

fce_travel_default <- tribble(
  ~species, ~travel_form,
  "Chinook", "log-normal",
  "Coho", "exponential",
  "Steelhead", "exponential"
) |>
  mutate(
    travel_prior = list(c(10, 1), 4, 4)
  )

fce_gen_pcap_prior <- function(n, uc_pcap_basis, sds, simplify = TRUE, ...) {
  nb <- ncol(uc_pcap_basis)
  length(sds) == nb || stop("Need one sd for each basis")
  replicate(n, uc_pcap_basis %*% rnorm(nb, 0, sds),
    simplify = simplify, ...
  )
}

## Full model functions --------------------------------------------------------
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Construct the data argument for the full model
##' @param pcap_formula A \emph{one-sided} formula potentially including both
##'   linear terms and mgcv-style smooth terms
##' @param pcap_sigma A named vector or named list with entries corresponding to
##'   the prior scale parameters for each linear term or smooth. Smooths should
##'   also include a \code{_null} term if any of the components (besides the
##'   last) are in the null space of the penalty.
##' @param rs_formula A \emph{one-sided} formula potentially including both
##'   linear terms and mgcv-style smooth terms
##' @param rs_sigma A named vector or named list with entries corresponding to
##'   the prior scale parameters for each linear term or smooth. Smooths should
##'   also include a \code{_null} term if any of the components (besides the
##'   last) are in the null space of the penalty.
##' @param release_df Data frame with one row per release
##' @param recap_df Data frame with one row per recapture opportunity, including
##'   all covariates used in \code{pcap_formula}
##' @param unmarked Data frame with one row per unmarked capture day, including
##'   all covariates used in \code{rs_formula}
##' @param trap_op Not currently used
##' @param travel_form Data frame indicating the form of the travel time
##'   distribution for each species. Distributions that can be used include the
##'   normal, inverse Gaussian, or exponential. Accepts multiple variations of
##'   these names.
##' @param travel_prior Expected number of days to arrive at the trap for each
##'   species. Order determined by species indexing (usually alphabetical).
##' @param validate Logical. Check dimensions etc? Turning this off can help
##'   with debugging.
##' @return A named list suitable for use in the FCE model.
##' @author John Best
##' @export
fce_data <- function(
    avail_formula, avail_sigma,
    pcap_formula, pcap_sigma,
    rs_formula, rs_sigma,
    release_df, recap_df, unmarked, trap_op,
    knots = NULL,
    validate = TRUE) {
  ## Join recapture observations with trap operations information
  recap_df <- recap_df |>
    left_join(select(trap_op, date, op),
      by = join_by(recapture_date == date)
    )

  ## Extract the start index of each release
  recap_idx <- recap_df |>
    summarize(
      n_recap = n(),
      .by = c(species, lifestage, release_date)
    )

  release_df <- left_join(
    release_df,
    recap_idx,
    by = join_by(species, lifestage, release_date)
  )

  recap_idx <- fce_ragged_index(release_df$n_recap)

  avail_pform <- fce_parse_formula(avail_formula, recap_df)
  avail_basis <- fce_construct_mm(avail_pform, recap_df, check_fullrank = FALSE)

  ## Construct the coefficient scaling parameter vector
  if (!all(unique(attr(avail_basis, "scale_index") %in% names(avail_sigma)))) {
    miss <- setdiff(attr(avail_basis, "scale_index"), names(avail_sigma))
    stop("avail_sigma missing value for ", miss)
  }
  avail_scale <- unlist(avail_sigma)[attr(avail_basis, "scale_index")]

  ## Construct smoother and design matrix for probability of capture
  pcap_pform <- fce_parse_formula(pcap_formula, recap_df, knots)
  pcap_basis <- fce_construct_mm(pcap_pform,
    check_fullrank = FALSE
  )
  rs_pcap_pred <- fce_construct_mm(pcap_pform,
    unmarked,
    check_fullrank = FALSE
  )

  ## Construct the coefficient scaling parameter vector
  if (!all(unique(attr(pcap_basis, "scale_index") %in% names(pcap_sigma)))) {
    miss <- setdiff(attr(pcap_basis, "scale_index"), names(pcap_sigma))
    stop("pcap_sigma missing value for ", miss)
  }
  pcap_scale <- unlist(pcap_sigma)[attr(pcap_basis, "scale_index")]

  ## Construct smoother and design matrix for runsize estimation
  rs_pform <- fce_parse_formula(rs_formula, unmarked)
  rs_basis <- fce_construct_mm(rs_pform,
    check_fullrank = FALSE
  )

  if (!all(unique(attr(rs_basis, "scale_index") %in% names(rs_sigma)))) {
    miss <- setdiff(attr(rs_basis, "scale_index"), names(rs_sigma))
    stop("rs_sigma missing value for ", miss)
  }
  rs_scale <- unlist(rs_sigma)[attr(rs_basis, "scale_index")]

  dat <- list(
    ## Sizes
    N_release = nrow(release_df),
    N_avail_basis = ncol(avail_basis),
    N_pcap_basis = ncol(pcap_basis),
    N_rec_obs = nrow(recap_df),
    N_runsize_basis = ncol(rs_basis),
    N_unmarked_obs = nrow(unmarked),

    ## Indices
    rec_idx = as.array(recap_idx),

    ## Release information
    num_released = release_df$count,
    num_lost = fce_num_lost(release_df, recap_df),

    ## Availability model
    avail_basis = avail_basis,
    avail_sigma = avail_scale,

    ## Recapture information
    num_recaptured = recap_df$count,
    rec_uc_pcap_basis = pcap_basis,
    pcap_sigma = pcap_scale,

    ## Unmarked capture observations
    unm_count = unmarked$count,
    runsize_basis = rs_basis,
    runsize_sigma = rs_scale,
    runsize_pcap_pred = rs_pcap_pred
  )

  ## Save the formula specifications so they can be used for plotting later
  attr(dat, "avail_pform") <- avail_pform
  attr(dat, "pcap_pform") <- pcap_pform
  attr(dat, "rs_pform") <- rs_pform
  attr(dat, "avail_scale_index") <- attr(avail_basis, "scale_index")
  attr(dat, "pcap_scale_index") <- attr(pcap_basis, "scale_index")
  attr(dat, "rs_scale_index") <- attr(rs_basis, "scale_index")

  ## Helpful to turn off validation to track down errors (maybe this should be a
  ## warning instead?)
  if (validate) {
    fce_validate_data(dat)
  }

  return(dat)
}

fce_validate_data <- function(data) {
  stopifnot(
    ## Sizes
    data$N_release > 0,
    data$N_avail_basis > 0,
    data$N_pcap_basis > 0,
    data$N_rec_obs > 0,
    data$N_runsize_basis > 0,
    data$N_unmarked_obs > 0,

    ## Indices
    length(data$rec_idx) == data$N_release + 1,
    all(data$rec_idx > 0),
    !is.unsorted(data$rec_idx),

    ## Release information
    length(data$num_released) == data$N_release,
    length(data$num_lost) == data$N_release,
    all(data$num_lost <= data$num_released),

    ## Availability model
    nrow(data$avail_basis) == data$N_rec_obs,
    ncol(data$avail_basis) == data$N_avail_basis,
    length(data$avail_sigma) == data$N_avail_basis,
    all(data$avail_sigma > 0),

    ## Recapture information
    length(data$num_recaptured) == data$N_rec_obs,
    nrow(data$rec_uc_pcap_basis) == data$N_rec_obs,
    ncol(data$rec_uc_pcap_basis) == data$N_pcap_basis,
    length(data$pcap_sigma) == data$N_pcap_basis,
    all(data$pcap_sigma > 0),

    ## Unmarked capture observations
    length(data$unm_count) == data$N_unmarked_obs,
    all(data$unm_count >= 0),
    nrow(data$runsize_basis) == data$N_unmarked_obs,
    ncol(data$runsize_basis) == data$N_runsize_basis,
    length(data$runsize_sigma) == data$N_runsize_basis,
    all(data$runsize_sigma > 0),
    nrow(data$runsize_pcap_pred) == data$N_unmarked_obs,
    ncol(data$runsize_pcap_pred) == data$N_pcap_basis
  )
  invisible(data)
}

fce_partial_effect <- function(use_pars, newdata, post, data, mod = "pcap") {
  pform <- attr(data, paste0(mod, "_pform"))
  ## scl_idx <- attr(data, paste0(mod, "_scale_index"))
  ## Get coefficient posteriors and rescale them
  ## coef <- post[[paste0(mod, "_coef")]] * data[[paste0(mod, "_sigma")]]

  mm <- fce_construct_mm(pform, newdata, check_fullrank = FALSE)
  ## zero_cols <- which(!(colnames(mm) %in% use_pars))
  keep_cols <- which(colnames(mm) %in% use_pars)
  ## mm[, zero_cols] <- 0

  mod <- ifelse(mod == "rs", "runsize", mod)
  coef <- post[[paste0(mod, "_coef")]][keep_cols] * data[[paste0(mod, "_sigma")]][keep_cols]

  eff <- mm[, keep_cols] %**% coef
  eff[, 1, drop = TRUE]
}
