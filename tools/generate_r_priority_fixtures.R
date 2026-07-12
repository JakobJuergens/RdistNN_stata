# Regenerate the high-priority cross-language values used by
# tests/test-r-parity-priority.do. Run from the RdistNN_stata repository root.

r_root <- normalizePath(file.path("..", "RdistNN"), mustWork = TRUE)
devtools::load_all(r_root, quiet = TRUE)

fixture_lines <- character()

emit_matrix <- function(name, x) {
  x <- as.matrix(x)
  rows <- apply(x, 1L, function(row) {
    values <- ifelse(is.na(row), ".", sprintf("%.17g", row))
    paste(values, collapse = ", ")
  })
  fixture_lines <<- c(
    fixture_lines,
    sprintf("matrix R_%s = (%s)", name, paste(rows, collapse = " \\ "))
  )
}

X <- matrix(
  c(
    -1, 0,
    0, -1,
    0, 0,
    1, 0,
    0, 1,
    1, 1,
    -1, -1,
    2, -0.5,
    -0.5, 1.5,
    0.75, -1.25,
    1.5, 0.5,
    -1.5, 0.75
  ),
  ncol = 2L,
  byrow = TRUE,
  dimnames = list(NULL, c("x1", "x2"))
)
n <- nrow(X)
i <- seq_len(n)
x0 <- c(0.2, -0.1)
X0 <- rbind(x0, c(-0.5, 0.75), c(1.2, 0.4))
y_np <- 0.7 * X[, 1L] - 0.4 * X[, 2L] + i / 20
rho <- 1 + 0.3 * X[, 1L] - 0.2 * X[, 2L] + (-1)^i * 0.25

D <- as.integer(i %% 3L != 0L)
e <- 0.2 + 0.05 * (i - 1L)
m0 <- 0.5 + 0.2 * X[, 1L] - 0.1 * X[, 2L]
m1 <- m0 + 1 + 0.1 * X[, 1L]
y_cate <- m0 + D * (m1 - m0) + ((i %% 4L) - 1.5) / 10
m <- e * m1 + (1 - e) * m0

md1 <- 0.4 * X[, 1L] - 0.2 * X[, 2L]
md2 <- -0.1 * X[, 1L] + 0.3 * X[, 2L]
d1 <- md1 + i / 30
d2 <- md2 - i / 40
g <- 0.15 + 0.1 * X[, 1L]
theta <- c(0.6, -0.4)
y_plqr <- d1 * theta[[1L]] + d2 * theta[[2L]] + g + ((i %% 5L) - 2) / 10

dat <- data.frame(
  x1 = X[, 1L], x2 = X[, 2L],
  y_np = y_np, rho = rho,
  D = D, e = e, m0 = m0, m1 = m1, m = m, y_cate = y_cate,
  d1 = d1, d2 = d2, md1 = md1, md2 = md2, g = g, y_plqr = y_plqr
)

for (method in c("exact", "asymp")) {
  np_fit <- np_dnn(
    X = X, x0 = x0, Y = y_np,
    method = method, k = 2L, s = 5L,
    se = TRUE, jackknife = "fast", return = "ess"
  )
  emit_matrix(paste0("np_", method), c(np_fit$mu, np_fit$se, np_fit$ess))

  cate_fit <- cate_dnn(
    X = X, X0 = X0, pseudo_outcome = rho,
    method = method, k = 2L, s = 5L,
    se = TRUE, jackknife = "fast", return = c("W", "ess")
  )
  emit_matrix(paste0("cate_", method, "_tau"), cate_fit$tau)
  emit_matrix(paste0("cate_", method, "_se"), cate_fit$se)
  emit_matrix(paste0("cate_", method, "_ess"), cate_fit$ess)
  cate_paper_diag <- diagnostics(cate_fit)[, c(
    "physical_radius_90", "arm_ess_available", "treated_ess", "control_ess",
    "score_contribution_available", "top5_score_share"
  ), drop = FALSE]
  emit_matrix(
    paste0("cate_", method, "_paper_diagnostics"),
    data.matrix(cate_paper_diag)
  )

  plqr_fit <- plqr_dnn(
    X = X, D = cbind(d1, d2), g = g, mD = cbind(md1, md2),
    Y = y_plqr, theta_eval = theta, tau = 0.5, x0 = x0,
    method = method, k = 2L, s = 5L,
    se = TRUE, jackknife = "fast", return = "ess"
  )
  plqr_score <- .plqr_dnn_score_matrix(
    cbind(d1, d2), cbind(md1, md2), g, y_plqr,
    theta_eval = theta, tau = 0.5
  )
  plqr_loo <- .dnn_jackknife_loo_array(
    x_mat = X, x0_mat = matrix(x0, nrow = 1L), signal = plqr_score,
    distance = "l2", method = method, k = 2L, s = 5L,
    normalize = TRUE, ctx = NULL, signal_arg = "score", jackknife = "fast"
  )
  plqr_vcov <- .dnn_jackknife_vcov_from_loo(
    loo = plqr_loo, q = 2L, n = n, m = 1L
  )$vcov[[1L]]
  emit_matrix(paste0("plqr_", method, "_moment"), plqr_fit$moment)
  emit_matrix(paste0("plqr_", method, "_moment_vcov"), plqr_vcov)
  emit_matrix(paste0("plqr_", method, "_se"), plqr_fit$moment_se)
  emit_matrix(paste0("plqr_", method, "_ess"), plqr_fit$ess)
  emit_matrix(paste0("plqr_", method, "_theta_eval"), plqr_fit$theta_eval)
}

np_batch_fit <- np_dnn(
  X = X, X0 = X0, Y = y_np,
  method = "exact", k = 2L, s = 5L,
  se = TRUE, jackknife = "fast", return = c("W", "ess", "dx")
)
np_batch_diag <- diagnostics(np_batch_fit)
np_batch_diag_cols <- c(
  "n_eff", "max_weight", "top1_weight_share", "top5_weight_share",
  "top10_weight_share", "n_nonzero_weight",
  "distance_diagnostics_available", "radius_rank_k", "radius_rank_s",
  "weighted_mean_distance", "weighted_median_distance",
  "outcome_available", "weighted_mean_y", "weighted_sd_y",
  "weighted_y_q10", "weighted_y_q50", "weighted_y_q90"
)
emit_matrix("np_batch_exact_mu", np_batch_fit$mu)
emit_matrix("np_batch_exact_se", np_batch_fit$se)
emit_matrix("np_batch_exact_ess", np_batch_fit$ess)
emit_matrix(
  "np_batch_exact_diagnostics",
  data.matrix(np_batch_diag[, np_batch_diag_cols, drop = FALSE])
)

np_grid_fit <- np_dnn(
  X = X, X0 = X0, Y = y_np,
  method = "exact", k = 2L,
  tuning_mode = "grid", se = TRUE, jackknife = "fast",
  return = c("W", "ess", "dx")
)
np_grid <- np_grid_fit$tuning_grid
np_grid_candidates <- cbind(
  candidate_id = seq_len(nrow(np_grid$candidates)),
  computed_id = match(
    np_grid$candidates$computed_key,
    unique(np_grid$candidates$computed_key)
  ),
  is_default = as.integer(np_grid$candidates$is_default),
  feasible = as.integer(np_grid$candidates$feasible),
  reason = ifelse(np_grid$candidates$reason == "ok", 0, 1),
  status = as.integer(np_grid$candidates$status == "ok"),
  lambda = np_grid$candidates$lambda,
  s = np_grid$candidates$s
)
np_grid_est <- matrix(
  np_grid$estimates$estimate,
  nrow = nrow(np_grid$candidates), byrow = TRUE
)
np_grid_se <- matrix(
  np_grid$estimates$se,
  nrow = nrow(np_grid$candidates), byrow = TRUE
)
np_grid_diag <- data.matrix(np_grid$diagnostics[, c(
  "ess", "max_weight", "top_n", "top5_weight_share", "degenerate"
), drop = FALSE])
emit_matrix("np_grid_candidates", np_grid_candidates)
emit_matrix("np_grid_estimates", np_grid_est)
emit_matrix("np_grid_se", np_grid_se)
emit_matrix("np_grid_default_ess", np_grid_fit$ess)
emit_matrix("np_grid_diagnostics", np_grid_diag)

cate_grid_fit <- cate_dnn(
  X = X, X0 = X0, pseudo_outcome = rho,
  method = "exact", k = 2L,
  tuning_mode = "grid", se = TRUE, jackknife = "fast",
  return = "ess"
)
cate_grid <- cate_grid_fit$tuning_grid
cate_grid_candidates <- cbind(
  candidate_id = seq_len(nrow(cate_grid$candidates)),
  computed_id = match(
    cate_grid$candidates$computed_key,
    unique(cate_grid$candidates$computed_key)
  ),
  is_default = as.integer(cate_grid$candidates$is_default),
  feasible = as.integer(cate_grid$candidates$feasible),
  reason = ifelse(cate_grid$candidates$reason == "ok", 0, 1),
  status = as.integer(cate_grid$candidates$status == "ok"),
  lambda = cate_grid$candidates$lambda,
  s = cate_grid$candidates$s
)
cate_grid_est <- matrix(
  cate_grid$estimates$estimate,
  nrow = nrow(cate_grid$candidates), byrow = TRUE
)
cate_grid_se <- matrix(
  cate_grid$estimates$se,
  nrow = nrow(cate_grid$candidates), byrow = TRUE
)
cate_grid_diag <- data.matrix(cate_grid$diagnostics[, c(
  "ess", "max_weight", "top_n", "top5_weight_share", "degenerate"
), drop = FALSE])
emit_matrix("cate_grid_candidates", cate_grid_candidates)
emit_matrix("cate_grid_estimates", cate_grid_est)
emit_matrix("cate_grid_se", cate_grid_se)
emit_matrix("cate_grid_default_ess", cate_grid_fit$ess)
emit_matrix("cate_grid_diagnostics", cate_grid_diag)

aipw_fit <- cate_dnn(
  data = dat, x_cols = c("x1", "x2"), y_col = "y_cate", t_col = "D",
  m0_col = "m0", m1_col = "m1", e_col = "e", x0 = x0,
  method = "exact", k = 2L, s = 5L,
  se = TRUE, jackknife = "fast", return = c("W", "ess")
)
robinson_fit <- cate_dnn(
  data = dat, x_cols = c("x1", "x2"), y_col = "y_cate", t_col = "D",
  m_col = "m", e_col = "e", x0 = x0,
  method = "exact", k = 2L, s = 5L,
  se = TRUE, jackknife = "fast", return = "ess"
)
emit_matrix("cate_aipw_exact", c(aipw_fit$tau, aipw_fit$se, aipw_fit$ess))
emit_matrix(
  "cate_aipw_paper_diagnostics",
  data.matrix(diagnostics(aipw_fit)[, c(
    "physical_radius_90", "arm_ess_available", "treated_ess", "control_ess",
    "score_contribution_available", "top5_score_share"
  ), drop = FALSE])
)
emit_matrix(
  "cate_robinson_exact",
  c(robinson_fit$tau, robinson_fit$se, robinson_fit$ess)
)

output <- file.path("tests", "fixtures", "r_priority_fixtures.do")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
header <- c(
  "// Generated by tools/generate_r_priority_fixtures.R; do not edit manually.",
  sprintf(
    "// RdistNN version %s; deterministic shared-surface parity fixture.",
    as.character(utils::packageVersion("RdistNN"))
  ),
  "version 17",
  ""
)
writeLines(c(header, fixture_lines), output, useBytes = TRUE)
cat(sprintf("Wrote %s with %d matrices.\n", output, length(fixture_lines)))
