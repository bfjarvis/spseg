# Numerical checks supporting the discussion draft; no package dependencies.
# Run from this directory with Rscript check-identities.R.

entropy <- function(p) {
  ans <- numeric(length(p))
  inside <- p > 0 & p < 1
  ans[inside] <- -p[inside] * log(p[inside]) -
    (1 - p[inside]) * log1p(-p[inside])
  ans
}
phi <- list(information = entropy,
            diversity = function(p) p * (1 - p),
            root = function(p) sqrt(p * (1 - p)))
denominators <- c(information = 1 / 2, diversity = 1 / 6, root = pi / 8)
models <- list(
  uniform = list(G = function(q) q, g = function(q) rep(1, length(q))),
  increasing = list(G = function(q) q^2, g = function(q) 2 * q),
  decreasing = list(G = function(q) 2 * q - q^2, g = function(q) 2 - 2 * q),
  cubic = list(G = function(q) q^3, g = function(q) 3 * q^2),
  mixed = list(G = function(q) .3 * q + .7 * q^3,
               g = function(q) .3 + 2.1 * q^2)
)

# Evaluate the triangular pair integral by midpoint quadrature. Prefix sums
# evaluate all ordered pairs without allocating a quadratic-size matrix.
pair_integral <- function(model, variation, n) {
  q <- (seq_len(n) - .5) / n
  p <- model$G(q)
  a <- variation(p) / (p * (1 - p))
  primitive <- (cumsum(a) - a / 2) / n
  mass <- model$g(q) / n
  left_mass <- c(0, head(cumsum(mass), -1))
  left_primitive <- c(0, head(cumsum(mass * primitive), -1))
  sum(mass * (primitive * left_mass - left_primitive))
}

checks <- do.call(rbind, lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  do.call(rbind, lapply(names(phi), function(family) {
    f <- phi[[family]]
    threshold <- integrate(function(q) f(model$G(q)), 0, 1,
                           rel.tol = 1e-10)$value
    pair <- pair_integral(model, f, 20000)
    data.frame(model = model_name, family = family, threshold = threshold,
               pair = pair, error = abs(pair - threshold))
  }))
}))
stopifnot(max(checks$error) < 1e-6)
for (family in names(phi)) {
  numerical <- integrate(phi[[family]], 0, 1, rel.tol = 1e-10)$value
  stopifnot(abs(numerical - denominators[[family]]) < 1e-9)
}

# Independently integrate the two arrangements in the paper's table.
cdf_from_intervals <- function(q, starts, ends) {
  result <- numeric(length(q))
  for (i in seq_along(starts)) {
    result <- result + pmax(0, pmin(q, ends[i]) - starts[i])
  }
  result / sum(ends - starts)
}
arrangements <- list(
  contiguous = list(list(starts = 0, ends = .5),
                    list(starts = .5, ends = 1)),
  interleaved = list(list(starts = c(0, .5), ends = c(.25, .75)),
                     list(starts = c(.25, .75), ends = c(.5, 1)))
)
table_values <- do.call(rbind, lapply(names(arrangements), function(name) {
  values <- vapply(names(phi), function(family) {
    local <- vapply(arrangements[[name]], function(intervals) {
      sum(vapply(seq_len(4), function(b) {
        integrate(function(q) phi[[family]](cdf_from_intervals(
          q, intervals$starts, intervals$ends)), (b - 1) / 4, b / 4,
          rel.tol = 1e-10)$value
      }, numeric(1)))
    }, numeric(1))
    1 - mean(local) / denominators[[family]]
  }, numeric(1))
  data.frame(arrangement = name, t(values), slice_average = .5)
}))
expected <- c(.5 - log(2) / 2, .125, .5 - 1 / pi)
stopifnot(max(abs(as.numeric(table_values[1, 2:4]) - .5)) < 1e-9)
stopifnot(max(abs(as.numeric(table_values[2, 2:4]) - expected)) < 1e-9)

# Existing population mixture: .5*q^2 + .5*(2*q-q^2) = q.
# Reardon's relative diversity is then 6*integral((q^2-q)^2) = .2.
reardon <- 6 * integrate(function(q) (q^2 - q)^2, 0, 1)$value
pair_result <- 1 - mean(vapply(models[c("increasing", "decreasing")],
  pair_integral, numeric(1), variation = phi$diversity, n = 20000)) / (1 / 6)
stopifnot(abs(reardon - .2) < 1e-10, abs(pair_result - reardon) < 1e-6)

# Check the slice/pair change of measure in both directions, including
# neighborhoods with zero density on parts of the rank interval.
xlogx <- function(x) {
  ans <- numeric(length(x))
  positive <- x > 0
  ans[positive] <- x[positive] * log(x[positive])
  ans
}
mixing <- list(
  information = function(x, y) (xlogx(x + y) - xlogx(x) - xlogx(y)) /
    (2 * log(2)),
  diversity = function(x, y) {
    ans <- numeric(length(x))
    positive <- x + y > 0
    ans[positive] <- 2 * x[positive] * y[positive] /
      (x[positive] + y[positive])
    ans
  },
  root = function(x, y) sqrt(x * y)
)
interval_model <- function(intervals) {
  list(
    G = function(q) cdf_from_intervals(q, intervals$starts, intervals$ends),
    g = function(q) {
      ans <- numeric(length(q))
      for (i in seq_along(intervals$starts)) {
        ans <- ans + (q > intervals$starts[i] & q < intervals$ends[i])
      }
      ans / sum(intervals$ends - intervals$starts)
    }
  )
}
mixtures <- c(list(smooth = models[c("increasing", "decreasing")]),
              lapply(arrangements, function(x) lapply(x, interval_model)))
n_pair <- 600L
q_pair <- (seq_len(n_pair) - .5) / n_pair
ij <- which(outer(q_pair, q_pair, "<"), arr.ind = TRUE)
r <- q_pair[ij[, 1]]
s <- q_pair[ij[, 2]]
weights <- list(equal = rep(1, length(r)), gap = s - r,
                lower_ranks = 2 - r - s)
q_fine <- (seq_len(20000L) - .5) / 20000L

bridge_checks <- list()
for (family in names(phi)) {
  f <- phi[[family]]
  B <- mixing[[family]]
  # Independently compare the closed form with the defining perspective.
  x <- c(0, 0, .2, 1, 2)
  y <- c(0, 2, 1.4, 1, 0)
  definition <- numeric(length(x))
  positive <- x + y > 0
  definition[positive] <- (x[positive] + y[positive]) / 2 *
    f(x[positive] / (x[positive] + y[positive])) / f(.5)
  stopifnot(max(abs(B(x, y) - definition)) < 1e-12)

  # Compute Reardon's A independently from its threshold-level primitive.
  pair_components <- function(model) {
    p <- model$G(q_fine)
    a <- numeric(length(p))
    inside <- p > 0 & p < 1
    a[inside] <- f(p[inside]) / (p[inside] * (1 - p[inside]))
    primitive <- (cumsum(a) - a / 2) / length(q_fine)
    at_pairs <- approx(q_fine, primitive, xout = q_pair)$y
    list(A = at_pairs[ij[, 2]] - at_pairs[ij[, 1]],
         x = model$g(r), y = model$g(s))
  }
  overall <- pair_components(models$uniform)
  reference <- sum(overall$A)
  for (name in names(mixtures)) {
    local <- lapply(mixtures[[name]], pair_components)
    B_local <- lapply(local, function(z) B(z$x, z$y))
    L <- 1 - Reduce(`+`, B_local) / length(local)
    for (K in weights) {
      converted <- vapply(seq_along(local), function(j) {
        z <- local[[j]]
        positive <- z$x * z$y > 0
        A_slice <- numeric(length(r))
        A_slice[positive] <- K[positive] * B_local[[j]][positive] /
          (z$x[positive] * z$y[positive])
        sum(A_slice * z$x * z$y)
      }, numeric(1))
      stopifnot(abs(sum(K * L) / sum(K) -
                      (1 - mean(converted) / sum(K))) < 1e-12)
    }
    weighted_mixing <- vapply(seq_along(local), function(j) {
      z <- local[[j]]
      positive <- z$x * z$y > 0
      W <- numeric(length(r))
      W[positive] <- z$A[positive] * z$x[positive] * z$y[positive] /
        B_local[[j]][positive]
      sum(W * B_local[[j]])
    }, numeric(1))
    recovered <- 1 - mean(weighted_mixing) / reference
    threshold_local <- vapply(mixtures[[name]], function(model) {
      sum(vapply(seq_len(4), function(b) {
        integrate(function(q) f(model$G(q)), (b - 1) / 4, b / 4,
                  rel.tol = 1e-10)$value
      }, numeric(1)))
    }, numeric(1))
    expected <- 1 - mean(threshold_local) / denominators[[family]]
    bridge_checks[[length(bridge_checks) + 1L]] <- data.frame(
      mixture = name, family = family, threshold = expected,
      weighted_slices = recovered, error = abs(recovered - expected))
  }
}
bridge_checks <- do.call(rbind, bridge_checks)
stopifnot(max(bridge_checks$error) < 1e-4)

# Changing the income schedule can alter inequality while preserving ranks.
income_schedules <- list(uniform = function(r) r, cubic = function(r) r^3)
income_ginis <- vapply(income_schedules, function(Q) {
  numerator <- integrate(function(r) (2 * r - 1) * Q(r), 0, 1)$value
  denominator <- integrate(Q, 0, 1)$value
  stopifnot(identical(rank(Q(q_pair)), rank(q_pair)))
  numerator / denominator
}, numeric(1))
stopifnot(max(abs(income_ginis - c(1 / 3, 3 / 5))) < 1e-10)

# Beta-product normalization and the analytic quarter-separation example.
shapes <- list(U = c(1, 1), P = c(1, 4), R = c(4, 1), M = c(4, 4))
beta_pairs <- list(UU = c("U", "U"), PU = c("P", "U"),
                   RU = c("R", "U"), PR = c("P", "R"), MM = c("M", "M"))
beta_values <- do.call(rbind, lapply(names(beta_pairs), function(name) {
  shape_f <- shapes[[beta_pairs[[name]][1]]]
  shape_h <- shapes[[beta_pairs[[name]][2]]]
  f <- function(x) dbeta(x, shape_f[1], shape_f[2])
  h <- function(x) dbeta(x, shape_h[1], shape_h[2])
  F <- function(x) pbeta(x, shape_f[1], shape_f[2])
  H <- function(x) pbeta(x, shape_h[1], shape_h[2])
  triangular_mass <- integrate(function(r) {
    (f(r) * (1 - H(r)) + h(r) * (1 - F(r))) / 2
  }, 0, 1, rel.tol = 1e-10)$value
  stopifnot(abs(triangular_mass - .5) < 1e-10)
  values <- vapply(c(.25, .75), function(cut) {
    # Independently integrate the cross-neighborhood rectangle.
    upper_f <- integrate(f, cut, 1, rel.tol = 1e-10)$value
    upper_h <- integrate(h, cut, 1, rel.tol = 1e-10)$value
    numerical <- integrate(function(r) {
      (f(r) * upper_h + h(r) * upper_f) / 2
    }, 0, cut, rel.tol = 1e-10)$value / triangular_mass
    analytic <- F(cut) * (1 - H(cut)) + H(cut) * (1 - F(cut))
    stopifnot(abs(numerical - analytic) < 1e-10)
    analytic
  }, numeric(1))
  data.frame(weight = name, bottom_quarter = values[1], top_quarter = values[2])
}))
stopifnot(max(abs(beta_values$bottom_quarter -
  c(.375, .591796875, .251953125, .682159423828125, .1311568021774292))) < 1e-10)
stopifnot(beta_values$bottom_quarter[2] > beta_values$top_quarter[2],
          beta_values$bottom_quarter[3] < beta_values$top_quarter[3])
for (cut in c(.25, .75)) {
  for (family in names(phi)) {
    f <- phi[[family]]
    lower <- integrate(function(q) f(q / cut), 0, cut,
                       rel.tol = 1e-10)$value
    upper <- integrate(function(q) f((q - cut) / (1 - cut)), cut, 1,
                       rel.tol = 1e-10)$value
    index <- 1 - (cut * lower + (1 - cut) * upper) / denominators[[family]]
    stopifnot(abs(index - .375) < 1e-9)
  }
}

# Rank-cell decompositions use exact beta masses and constant slice values
# within each quartile. Compare with independent quadrature of the full surface.
cell_decomposition <- function(local, B, cuts, shape_f, shape_h) {
  mid <- (head(cuts, -1) + tail(cuts, -1)) / 2
  mass_f <- diff(pbeta(cuts, shape_f[1], shape_f[2]))
  mass_h <- diff(pbeta(cuts, shape_h[1], shape_h[2]))
  cells <- which(upper.tri(matrix(0, length(mid), length(mid)), diag = TRUE),
                 arr.ind = TRUE)
  a <- cells[, 1]
  b <- cells[, 2]
  omega <- mass_f[a] * mass_h[b] + mass_h[a] * mass_f[b]
  omega[a == b] <- omega[a == b] / 2
  D <- 1 - Reduce(`+`, lapply(local, function(model) {
    B(model$g(mid[a]), model$g(mid[b]))
  })) / length(local)
  stopifnot(abs(sum(omega) - 1) < 1e-12)
  data.frame(a, b, omega, D, C = omega * D)
}
decomposition_checks <- list()
for (family in names(mixing)) {
  B <- mixing[[family]]
  for (name in names(beta_pairs)) {
    shape_f <- shapes[[beta_pairs[[name]][1]]]
    shape_h <- shapes[[beta_pairs[[name]][2]]]
    K <- (dbeta(r, shape_f[1], shape_f[2]) * dbeta(s, shape_h[1], shape_h[2]) +
          dbeta(r, shape_h[1], shape_h[2]) * dbeta(s, shape_f[1], shape_f[2])) / 2
    levels <- list()
    for (arrangement in names(arrangements)) {
      local <- mixtures[[arrangement]]
      coarse <- cell_decomposition(local, B, seq(0, 1, .25), shape_f, shape_h)
      fine <- cell_decomposition(local, B, seq(0, 1, .125), shape_f, shape_h)
      # Refining reporting bands preserves each coarse contribution.
      for (i in seq_len(nrow(coarse))) {
        inside <- ceiling(fine$a / 2) == coarse$a[i] &
          ceiling(fine$b / 2) == coarse$b[i]
        stopifnot(abs(sum(fine$C[inside]) - coarse$C[i]) < 1e-12)
      }
      L <- 1 - Reduce(`+`, lapply(local, function(model) {
        B(model$g(r), model$g(s))
      })) / length(local)
      quadrature <- 2 * sum(K * L) / n_pair^2
      stopifnot(abs(quadrature - sum(coarse$C)) < 1e-5)
      levels[[arrangement]] <- coarse
    }
    before <- levels$contiguous
    after <- levels$interleaved
    changes <- after$C - before$C
    stopifnot(max(abs(changes - before$omega * (after$D - before$D))) < 1e-12)
    if (name == "UU") {
      expected <- matrix(0, 4, 4)
      expected[1, 2] <- expected[3, 4] <- 1 / 8
      expected[1, 3] <- expected[2, 4] <- -1 / 8
      stopifnot(abs(sum(before$C) - .5) < 1e-12,
                abs(sum(after$C) - .5) < 1e-12,
                max(abs(changes - expected[cbind(before$a, before$b)])) < 1e-12)
    }
    decomposition_checks[[length(decomposition_checks) + 1L]] <- data.frame(
      family, weight = name, before = sum(before$C), after = sum(after$C),
      change = sum(changes))
  }
}
decomposition_checks <- do.call(rbind, decomposition_checks)

print(checks, row.names = FALSE, digits = 8)
print(table_values, row.names = FALSE, digits = 8)
print(bridge_checks, row.names = FALSE, digits = 8)
print(beta_values, row.names = FALSE, digits = 8)
print(income_ginis, digits = 8)
print(decomposition_checks, row.names = FALSE, digits = 8)
cat("All identity, normalization, illustration, slice/pair, beta-product, and decomposition checks passed.\n")
