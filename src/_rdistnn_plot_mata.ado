version 17
// Numeric plotting adapters; estimator and weight policy stay in core modules.
capture mata: mata drop rdistnn_plot_target_diag()
capture mata: mata drop rdistnn_plot_candidate_diag()
capture mata: mata drop rdistnn_plot_est_default()
capture mata: mata drop rdistnn_plot_est_grid()
capture mata: mata drop rdistnn_plot_tuning()
capture mata: mata drop rdistnn_plot_overlap()
capture mata: mata drop rdistnn_plot_rank_weights()
mata:
real matrix rdistnn_plot_target_diag(real matrix diagnostics,
                                              real matrix targets)
{
    if (rows(diagnostics) != rows(targets)) {
        errprintf("rdistnn_plot_target_diagnostics: row counts must match.\n")
        exit(error(503))
    }
    return(((1::rows(diagnostics)), targets, diagnostics))
}

real matrix rdistnn_plot_candidate_diag(real matrix meta,
                                                 real matrix diagnostics,
                                                 real colvector idx)
{
    if (rows(meta) != rows(diagnostics)) {
        errprintf("rdistnn_plot_candidate_diagnostics: row counts must match.\n")
        exit(error(503))
    }
    return((meta[idx,], diagnostics[idx,]))
}

real matrix rdistnn_plot_est_default(real matrix targets,
                                            real colvector estimate,
                                            real colvector se,
                                            real scalar s,
                                            real scalar z)
{
    real scalar m
    m = rows(targets)
    if (rows(estimate) != m | rows(se) != m) {
        errprintf("rdistnn_plot_estimates_default: estimate and se must match targets.\n")
        exit(error(503))
    }
    return((targets, J(m, 1, 0), J(m, 1, s), J(m, 1, 1),
        (1::m), estimate, se, estimate :- z :* se, estimate :+ z :* se))
}

real matrix rdistnn_plot_est_grid(real matrix targets,
                                         real matrix meta,
                                         real matrix estimates,
                                         real matrix ses,
                                         real colvector idx,
                                         real scalar z)
{
    real scalar i, g, m
    real colvector b, se
    real matrix out

    m = rows(targets)
    out = J(0, cols(targets) + 8, .)
    for (i = 1; i <= rows(idx); i++) {
        g = idx[i]
        b = estimates[g,]'
        se = ses[g,]'
        out = out \ (targets, J(m, 1, g), J(m, 1, meta[g,8]),
            J(m, 1, meta[g,3]), (1::m), b, se,
            b :- z :* se, b :+ z :* se)
    }
    return(out)
}

real matrix rdistnn_plot_tuning(real matrix meta,
                                      real matrix estimates,
                                      real matrix ses,
                                      real colvector candidate_idx,
                                      real colvector target_idx,
                                      real scalar z)
{
    real scalar i, j, g, t, b, se
    real matrix out
    out = J(0, 9, .)
    for (i = 1; i <= rows(candidate_idx); i++) {
        g = candidate_idx[i]
        for (j = 1; j <= rows(target_idx); j++) {
            t = target_idx[j]
            b = estimates[g,t]
            se = ses[g,t]
            out = out \ (g, meta[g,8], meta[g,7], meta[g,3], t,
                b, se, b - z * se, b + z * se)
        }
    }
    return(out)
}

real matrix rdistnn_plot_overlap(real matrix diagnostics,
                                       real matrix targets,
                                       real rowvector cols_use)
{
    if (rows(diagnostics) != rows(targets)) {
        errprintf("rdistnn_plot_overlap_data: row counts must match.\n")
        exit(error(503))
    }
    return((targets, (1::rows(diagnostics)), diagnostics[,cols_use]))
}

real matrix rdistnn_plot_rank_weights(real matrix meta,
                                       real colvector idx,
                                       real scalar n,
                                       real scalar k,
                                       real scalar default_s,
                                       string scalar method,
                                       real scalar use_meta)
{
    real scalar i, g, s, is_default
    real colvector w
    real matrix out
    out = J(0, 6, .)
    for (i = 1; i <= rows(idx); i++) {
        g = idx[i]
        if (use_meta) {
            s = meta[g,8]
            is_default = meta[g,3]
        }
        else {
            s = default_s
            is_default = 1
        }
        w = rdistnn_weights_from_ranks((1::n), k, s, method)
        out = out \ ((1::n), w, runningsum(w), J(n,1,g),
            J(n,1,s), J(n,1,is_default))
    }
    return(out)
}

end
