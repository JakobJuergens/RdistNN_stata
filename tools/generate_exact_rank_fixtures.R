#!/usr/bin/env Rscript
# Run from the Stata checkout with the sibling R package available.
args <- commandArgs(trailingOnly = TRUE)
r_root <- if (length(args)) args[1L] else "../RdistNN"
pkgload::load_all(r_root, compile = FALSE, quiet = TRUE)
cases <- expand.grid(n = c(8L, 20L, 257L), k = c(1L, 3L), s = c(3L, 7L))
rows <- lapply(seq_len(nrow(cases)), function(i) {
  a <- cases[i, ]
  options(RdistNN.use_compiled = FALSE)
  w <- rank_weight_eval(rank_weight_dnn_exact(a$k, a$s), seq_len(a$n))
  options(RdistNN.use_compiled = TRUE)
  stopifnot(is.function(RdistNN:::.rdistnn_compiled_fun("rustat_exact_rank_table")),
    identical(w, rank_weight_eval(rank_weight_dnn_exact(a$k, a$s), seq_len(a$n))))
  data.frame(n = a$n, k = a$k, s = a$s, rank = seq_len(a$n),
    weight = sprintf("%.17g", w))
})
write.table(do.call(rbind, rows), "tests/fixtures/exact_rank_tables.csv",
  sep = ",", row.names = FALSE, quote = FALSE)
