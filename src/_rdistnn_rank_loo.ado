*! version 0.1.5  2026-10-25
version 17
// Arbitrary-rank-table deletions and DNN covariance adapters; TDNN reuses tables.
foreach fn in naive_table prefix plugin from_ranks from_table {
    capture mata: mata drop rdistnn_loo_`fn'()
}
capture mata: mata drop rdistnn_center_signal()
capture mata: mata drop rdistnn_jk_cov_from_ranks()
capture mata: mata drop rdistnn_jk_cov_multi_from_ranks()
mata:
// Explicit table-index deletions are the reference and final fallback.
real matrix rdistnn_loo_naive_table(real matrix ranks, real matrix signal,
    real colvector table)
{
    real scalar n,m,q,i,j
    real colvector kept,r,active,ord
    real rowvector level
    real matrix loo,z
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    loo=J(n,m*q,.)
    active=((table:!=0)\0):|(0\(table:!=0))
    for(j=1;j<=m;j++) {
        ord=J(n,1,.)
        ord[ranks[,j]]=(1::n)
        z=rdistnn_center_signal(signal[ord,],active,level,table)
        z=z[ranks[,j],]
        for(i=1;i<=n;i++) {
            kept=((1::n):!=i)
            r=select(ranks[,j],kept)
            r=r:-(r:>ranks[i,j])
            loo[i,((j-1)*q+1)..(j*q)]=rdistnn_weighted_signals(table[r],select(z,kept)):+level*rdistnn_sum(table)
        }
    }
    return(loo)
}

// Two linear passes per target/component. The table at n-1 is shared across
// every deletion and target; no deletion-specific rank or signal copies.
real matrix rdistnn_loo_prefix(real matrix ranks, real matrix signal,
    real colvector table)
{
    real scalar n, m, q, j, mass
    real colvector ord, active
    real rowvector level
    real matrix z, values, loo
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    loo=J(n,m*q,.)
    mass=rdistnn_sum(table)
    for(j=1;j<=m;j++) {
        ord=J(n,1,.)
        ord[ranks[,j]]=(1::n)
        active=((table:!=0)\0):|(0\(table:!=0))
        z=rdistnn_center_signal(signal[ord,],active,level,table)
        values=rdistnn_accumulation_deletions(z,table,level*mass)
        loo[,((j-1)*q+1)..(j*q)]=values[ranks[,j],]
    }
    return(loo)
}

// The bridge uses Mata-created Stata matrices, which may exceed c(matsize).
// Attempt the native engine without changing the user's matrix setting.
real matrix rdistnn_loo_plugin(real matrix ranks, real matrix signal,
    real colvector table)
{
    if(st_global("rdistnn_plugin_ok")!="1" | !rdistnn_native_packed_api()) return(J(0,0,.))
    return(rdistnn_native_loo(ranks,signal,table))
}

// Accompanying DX enables exact-distance groups; supplied ranks are never
// recomputed. Without DX, deletions use the strict rank table.
// Common dispatcher: signed and nonnegative tables use identical validation,
// centering and native -> Mata -> explicit-deletion fallback semantics.
real matrix rdistnn_loo_from_table(real matrix ranks, real matrix signal,
    real colvector table, | real scalar center, string scalar engine, real matrix DX, real matrix prepared_groups)
{
    real scalar n,m,q,i,j,r
    real colvector seen,active,ord
    real matrix loo,z
    real rowvector level
    string scalar plugin_mode
    pointer(function) scalar prefix
    if(args()<7) prepared_groups=J(0,0,.)
    if(args()<4) center=0
    if(args()<5 | engine=="") engine=st_global("rdistnn_jackknife_engine")
    if(engine=="") engine="auto"
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    if(n<2 | m<1 | rows(signal)!=n | q<1 | hasmissing(signal) |
        hasmissing(ranks) | min(ranks)<1 | max(ranks)>n |
        any(ranks:!=floor(ranks)) | rows(table)!=n-1 |
        cols(table)!=1 | hasmissing(table) | (center!=0 & center!=1)) {
        rdistnn_fail(198,"Rank-table deletions require finite aligned signals, rank permutations and n-1 weights.")
    }
    for(j=1;j<=m;j++) {
        seen=J(n,1,0)
        for(i=1;i<=n;i++) {
            r=ranks[i,j]
            if(seen[r]) {
                rdistnn_fail(198,"Each ranks column must be a permutation of 1..n.")
            }
            seen[r]=1
        }
    }
    if(!any(engine:==("auto","plugin","mata","naive"))) {
        rdistnn_fail(198,"rdistnn_jackknife_engine must be auto, plugin, mata or naive.")
    }
    if(args()>=6 & rows(DX)>0) return(rdistnn_distance_loo(ranks,signal,table,center,engine,DX,prepared_groups))
    plugin_mode=st_global("rdistnn_use_plugin")
    if(engine=="auto" & plugin_mode!="1" & all(signal:==0)) {
        st_global("rdistnn_last_jackknife_engine","zero")
        return(J(n,m*q,0))
    }
    if(center) {
        active=((table:!=0)\0):|(0\(table:!=0))
        loo=J(n,m*q,.)
        for(j=1;j<=m;j++) {
            ord=J(n,1,.)
            ord[ranks[,j]]=(1::n)
            z=rdistnn_center_signal(signal[ord,],active,level,table)
            z=z[ranks[,j],]
            loo[,((j-1)*q+1)..(j*q)]=rdistnn_loo_from_table(ranks[,j],z,table,0,engine)
        }
        return(loo)
    }
    if(engine=="naive") {
        st_global("rdistnn_last_jackknife_engine","naive")
        return(rdistnn_loo_naive_table(ranks,signal,table))
    }
    if(engine=="plugin" | (engine=="auto" & plugin_mode!="0")) {
        loo=rdistnn_loo_plugin(ranks,signal,table)
        if(rows(loo)==n & cols(loo)==m*q & !hasmissing(loo)) {
            st_global("rdistnn_last_jackknife_engine","plugin")
            return(loo)
        }
        if(engine=="plugin" | plugin_mode=="1") {
            rdistnn_fail(498,"Compiled rank deletion API is unavailable; install an updated plugin or allow Mata fallback.")
        }
    }
    prefix=findexternal("rdistnn_loo_prefix()")
    if(prefix!=NULL) {
        loo=(*prefix)(ranks,signal,table)
        if(rows(loo)==n & cols(loo)==m*q & !hasmissing(loo)) {
            st_global("rdistnn_last_jackknife_engine","mata")
            return(loo)
        }
    }
    if(engine=="mata") {
        rdistnn_fail(498,"Mata prefix/suffix engine unavailable or produced nonfinite deletions; rescale signals or select another jackknife engine.")
    }
    st_global("rdistnn_last_jackknife_engine","naive")
    return(rdistnn_loo_naive_table(ranks,signal,table))
}

// DNN adapter constructs a deletion table; the dispatcher owns computation.
real matrix rdistnn_loo_from_ranks(real matrix ranks, real matrix signal,
    real scalar k, real scalar s, | string scalar method, real scalar center,
    real colvector precomputed_table, real matrix DX, real matrix prepared_groups)
{
    real scalar n
    real colvector table
    if(args()<9) prepared_groups=J(0,0,.)
    if(args()<5 | method=="") method="exact"
    if(args()<6) center=0
    n=rows(ranks)
    if((method!="exact" & method!="asymp") | n<2 | missing(k) |
        missing(s) | k<1 | k!=floor(k) | s<k | s>=n | s!=floor(s)) {
        rdistnn_fail(198,"rdistnn_loo_from_ranks: require 1 <= k <= s < n and exact/asymp weights.")
    }
    if(args()<7) table=rdistnn_weights_from_ranks((1::(n-1)),k,s,method)
    else table=precomputed_table
    if(args()>=8 & rows(DX)>0) return(rdistnn_loo_from_table(ranks,signal,table,center,"",DX,prepared_groups))
    return(rdistnn_loo_from_table(ranks,signal,table,center))
}

// A common signal level has zero jackknife covariance. Remove it before
// accumulating deletions, and never restore it on covariance-only paths.
// Use only the union of deletion support to choose the reference, and zero
// unused rows before subtraction. A deleted earlier rank can bring the next
// original rank into support. Signed tables require nonzero support.
real matrix rdistnn_center_signal(real matrix signal, | real colvector active,
    real rowvector level, real matrix weights)
{
    real matrix z, used
    if(args()<2) active=J(rows(signal),1,1)
    z=J(rows(signal),cols(signal),0)
    level=J(1,cols(signal),0)
    if(!any(active)) return(z)
    used=select(signal,active)
    if(args()>=4) level=rdistnn_signal_level(used,weights)
    else level=rdistnn_signal_level(used)
    z[selectindex(active),]=used:-level
    return(z)
}

// ── rdistnn_jk_cov_from_ranks ─────────────────────────────────────────
// Rank-shift jackknife covariance for a DNN-weighted signal matrix.  When
// observation i
// (with rank r_i) is deleted, all remaining observations with rank > r_i
// are shifted down by 1, preserving the relative ordering under a fixed
// distance metric.
//
// This is exact delete-1 rank recomputation for L2.  For Mahalanobis, callers
// interpret the covariance metric as fixed at the full estimation sample.
// This matches the "fast" jackknife in the R reference (jackknife_rank_full).
//
// Arguments:
//   ranks   : n-vector of integer ranks from the full sample (1 = nearest)
//   signal  : n x q matrix of per-observation contributions
//   k       : number selected within each s-subset
//   s       : kernel order
//   method  : "exact" (default) or "asymp"
//
// Returns: q x q jackknife covariance matrix.
//
real matrix rdistnn_jk_cov_from_ranks(real colvector ranks,
                                       real matrix signal,
                                       real scalar k,
                                       real scalar s,
                                       | string scalar method, real matrix DX, real matrix prepared_groups)
{
    real scalar     n, q
    real matrix     loo

    if(args()<7) prepared_groups=J(0,0,.)
    if (args() < 5 | method == "") method = "exact"

    n = rows(ranks)
    if (rows(signal) != n) {
        rdistnn_fail(503,"rdistnn_jk_cov_from_ranks: rows(signal) must equal rows(ranks).")
    }

    q = cols(signal)
    if (q < 1) {
        rdistnn_fail(198,"rdistnn_jk_cov_from_ranks: signal must have at least one column.")
    }

    if(args()>=6 & rows(DX)>0) {
        loo=rdistnn_loo_from_ranks(ranks,signal,k,s,method,1,
            rdistnn_weights_from_ranks((1::(n-1)),k,s,method),DX,prepared_groups)
    }
    else loo=rdistnn_loo_from_ranks(ranks,signal,k,s,method,1)
    if(all(loo:==0)) return(J(q,q,0))

    return(rdistnn_covariance(loo,(n-1)/n,1))
}

// Joint rank-shift jackknife covariance for a scalar signal evaluated at
// multiple targets. Each ranks column is the full-sample rank permutation for
// one target. The returned m x m matrix retains cross-target covariance.
real matrix rdistnn_jk_cov_multi_from_ranks(real matrix ranks,
                                             real colvector signal,
                                             real scalar k,
                                             real scalar s,
                                             | string scalar method, real matrix DX, real matrix prepared_groups)
{
    real scalar     n, m
    real matrix     loo

    if(args()<7) prepared_groups=J(0,0,.)
    if (args() < 5 | method == "") method = "exact"
    n = rows(ranks)
    m = cols(ranks)
    if (rows(signal) != n | cols(signal) != 1) {
        rdistnn_fail(503,"rdistnn_jk_cov_multi_from_ranks: signal must be an n x 1 vector matching ranks.")
    }
    if (m < 1) {
        rdistnn_fail(198,"rdistnn_jk_cov_multi_from_ranks: ranks must have at least one target column.")
    }

    if(args()>=6 & rows(DX)>0) {
        loo=rdistnn_loo_from_ranks(ranks,signal,k,s,method,1,
            rdistnn_weights_from_ranks((1::(n-1)),k,s,method),DX,prepared_groups)
    }
    else loo=rdistnn_loo_from_ranks(ranks,signal,k,s,method,1)
    if(all(loo:==0)) return(J(m,m,0))

    return(rdistnn_covariance(loo,(n-1)/n,1))
}

end
