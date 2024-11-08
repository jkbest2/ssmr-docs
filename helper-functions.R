### Helper functions for finding PC prior parameters etc.

scale_vec <- function(x, center = TRUE, scale = TRUE) {
  x_scl <- scale(x, center, scale)
  structure(c(x_scl),
    "scaled:center" = attr(x_scl, "scaled:center"),
    "scaled:scale" = attr(x_scl, "scaled:scale")
  )
}

scale_to <- function(x, scale) {
  x0 <- scale(x)
  xs <- scale * x0
  newscale <- attr(x0, "scaled:scale") / scale

  structure(
    c(xs),
    "scaled:center" = attr(x0, "scaled:center"),
    "scaled:scale" = newscale
  )
}

unscale_vec <- function(x_scl, center = NULL, scale = NULL) {
  if (is.null(center)) {
    center <- attr(x_scl, "scaled:center")
  }
  if (is.null(center)) {
    center <- 0
  }
  if (is.null(scale)) {
    scale <- attr(x_scl, "scaled:scale")
  }
  if (is.null(scale)) {
    scale <- 1
  }

  x_scl * scale + center
}

## Marginal precision PC prior parameter; Pr(marg std dev < U) = alpha
tau_lambda <- function(alpha, U) {
  alpha > 0 && alpha < 1 || stop("alpha must be in (0, 1)")
  -log(alpha) / U
}

## PC prior parameter to shrink a stationary AR(1) process towards constant;
## Pr(phi > U) = alpha. For theta = 10, U = 0.9, there is a roughly 95%
## probability that phi > 0.9
phi_theta <- function(alpha, U, uniroot_interval = c(0, 1000)) {
  alpha > 0 && alpha < 1 || stop("alpha must be in (0, 1)")
  alpha_fun <- function(theta) {
    1 - exp(-theta * sqrt(1 - U)) / (1 - exp(-sqrt(2) * theta)) - alpha
  }
  uniroot(alpha_fun, interval = uniroot_interval)$root
}

find_knots <- function(x, num_knots = 20, bdeg = 3, pad = 0.1) {
  xrange <- range(x) + c(-pad, pad) * sd(x)
  dx <- diff(xrange) / num_knots
  ## Extend the knot placement based on the number of knots, provides the extra
  ## knots depending on the degree of the spline
  seq(xrange[1] - bdeg * dx, xrange[2] + bdeg * dx, by = dx)
}

get_uc_pcap_pen_idx <- function(s) {
  idx <- diag(s)
  idx[idx == 0] <- c(2, 3)[seq_along(idx[idx == 0])]
  idx
}

str_extract_digits <- function(n) as.numeric(stringr::str_extract(n, "[:digit:]+"))

## Inverse Gaussian distribution
dinvgauss <- function(x, mu, lambda, log = FALSE) {
  stopifnot(x >= 0, mu > 0, lambda > 0)
  # if (x == 0) {
  #   ld <- -Inf
  # } else {
  ld1 <- log(lambda) / 2 - log(2 * pi) / 2 - 3 / 2 * log(x)
  ld2 <- -(lambda * (x - mu)^2) / (2 * mu^2 * x)
  ld <- ld1 + ld2
  # }
  if (!log) {
    ld <- exp(ld)
  }
  return(ld)
}

trajectory <- function(coef, basis, coef_scale = 1) {
  coef <- coef * coef_scale
  basis %*% coef
}

softmax <- function(x) {
  exp(x) / sum(exp(x))
}
