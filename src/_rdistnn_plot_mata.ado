*! version 0.1.5  2026-10-25
version 17
// Numeric plotting adapters; estimator and weight policy stay in core modules.
capture mata: mata drop rdistnn_plot_est_default()
capture mata: mata drop rdistnn_plot_est_grid()
capture mata: mata drop rdistnn_plot_tuning()
capture mata: mata drop rdistnn_plot_overlap()
capture mata: mata drop rdistnn_plot_rank_weights()
mata:



real matrix rdistnn_plot_est_default(real matrix targets,
                                            real colvector estimate,
                                            real colvector se,
                                            real scalar s,
                                            real scalar z)
{
    real scalar m
    m = rows(targets)
    if (rows(estimate) != m | rows(se) != m) {
        rdistnn_fail(503,"rdistnn_plot_estimates_default: estimate and se must match targets.")
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
        rdistnn_fail(503,"rdistnn_plot_overlap_data: row counts must match.")
    }
    return((targets, (1::rows(diagnostics)), diagnostics[,cols_use]))
}

real matrix rdistnn_plot_rank_weights(real matrix meta,
                                       real colvector idx,
                                       real scalar n,
                                       real scalar k,
                                       real scalar default_s,
                                       string scalar method,
                                       real scalar use_meta, real matrix ranks, real rowvector targets, real matrix DX)
{
    real scalar i, g, s, is_default,j,t
    real colvector w,ord
    real matrix out,W,groups,R,D
    out = J(rows(idx)*cols(targets)*n,7,.)
    R=ranks[,targets]; D=DX[,targets]
    groups=st_matrix("e(model_groups)")[,targets]
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
        W=rdistnn_weights_from_rmat(R,k,s,method,D,groups)
        for(j=1;j<=cols(targets);j++) {
            t=targets[j]; ord=J(n,1,.); ord[R[,j]]=(1::n)
            w=W[ord,j]
            out[((i-1)*cols(targets)+j-1)*n:+(1::n),]=((1::n),w,runningsum(w),J(n,1,g),
                J(n,1,s),J(n,1,is_default),J(n,1,t))
        }
    }
    return(out)
}

end

capture mata: mata drop rdistnn_plot_target_ids()
capture mata: mata drop rdistnn_plot_wald_available()
mata:
void rdistnn_plot_target_ids(string scalar data, real scalar column,
    string scalar destination)
{
    real colvector ids
    ids=st_matrix(data)[,column]
    st_local(destination,invtokens(strofreal(uniqrows(select(ids,ids:<.)))'))
}

real colvector rdistnn_plot_wald_available(real matrix work, real matrix status,
    real scalar target, real scalar candidate)
{
    real colvector available
    real scalar i,g,t,valid
    if(candidate==0) return(status[work[,target],5])
    available=J(rows(work),1,0)
    for(i=1;i<=rows(status);i++) {
        g=status[i,1]; t=status[i,2]; valid=status[i,7]==1
        available=available:+(work[,candidate]:==g):*(work[,target]:==t):*valid
    }
    return(available)
}
end
