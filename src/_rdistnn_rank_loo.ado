*! version 0.1.0  2026-10-18
version 17
// One arbitrary-rank-table deletion engine for DNN and signed TDNN adapters.
foreach fn in naive naive_table prefix plugin from_ranks from_table {
    capture mata: mata drop rdistnn_loo_`fn'()
}
capture mata: mata drop rdistnn_center_signal()
mata:
// Explicit table-index deletions are the reference and final fallback.
real matrix rdistnn_loo_naive_table(real matrix ranks, real matrix signal,
    real colvector table)
{
    real scalar n,m,q,i,j
    real colvector kept,r,active
    real rowvector level
    real matrix loo,z
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    loo=J(n,m*q,.)
    active=((table:!=0)\0):|(0\(table:!=0))
    for(j=1;j<=m;j++) {
        z=rdistnn_center_signal(signal,active[ranks[,j]],level)
        for(i=1;i<=n;i++) {
            kept=((1::n):!=i)
            r=select(ranks[,j],kept)
            r=r:-(r:>ranks[i,j])
            loo[i,((j-1)*q+1)..(j*q)]=cross(table[r],select(z,kept)):+level*sum(table)
        }
    }
    return(loo)
}

// Compatibility adapter for the original DNN reference API.
real matrix rdistnn_loo_naive(real matrix ranks, real matrix signal,
    real scalar k, real scalar s, string scalar method)
{
    return(rdistnn_loo_naive_table(ranks,signal,
        rdistnn_weights_from_ranks((1::(rows(ranks)-1)),k,s,method)))
}

// Two linear passes per target/component. The table at n-1 is shared across
// every deletion and target; no deletion-specific rank or signal copies.
real matrix rdistnn_loo_prefix(real matrix ranks, real matrix signal,
    real colvector table)
{
    real scalar n, m, q, j, r
    real colvector ord, active
    real rowvector level
    real matrix z, left, right, loo
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    loo=J(n,m*q,.)
    for(j=1;j<=m;j++) {
        ord=J(n,1,.)
        ord[ranks[,j]]=(1::n)
        active=((table:!=0)\0):|(0\(table:!=0))
        z=rdistnn_center_signal(signal[ord,],active,level)
        left=right=J(n,q,0)
        for(r=2;r<=n;r++) left[r,]=left[r-1,]+table[r-1]*z[r-1,]
        for(r=n-1;r>=1;r--) right[r,]=right[r+1,]+table[r]*z[r+1,]
        loo[,((j-1)*q+1)..(j*q)]=(left[ranks[,j],]+right[ranks[,j],]):+level*sum(table)
    }
    return(loo)
}

// The bridge uses Mata-created Stata matrices, which may exceed c(matsize).
// Attempt the native engine without changing the user's matrix setting.
real matrix rdistnn_loo_plugin(real matrix ranks, real matrix signal,
    real colvector table)
{
    string rowvector names
    real scalar rc, i
    real matrix loo
    if(st_global("rdistnn_plugin_ok")!="1") return(J(0,0,.))
    names=st_tempname(4)
    st_matrix(names[1],ranks)
    st_matrix(names[2],signal)
    st_matrix(names[3],table)
    st_matrix(names[4],J(rows(ranks),cols(ranks)*cols(signal),.))
    rc=_stata("plugin call rdistnn_plugin, rank_loo "+invtokens(names),1)
    loo=J(0,0,.)
    if(rc==0) loo=st_matrix(names[4])
    for(i=1;i<=4;i++) st_matrix(names[i],J(0,0,.))
    return(loo)
}

// Common dispatcher: signed and nonnegative tables use identical validation,
// centering and native -> Mata -> explicit-deletion fallback semantics.
real matrix rdistnn_loo_from_table(real matrix ranks, real matrix signal,
    real colvector table, | real scalar center, string scalar engine)
{
    real scalar n,m,q,i,j,r
    real colvector seen,active
    real matrix loo,z
    string scalar plugin_mode
    pointer(function) scalar prefix
    if(args()<4) center=0
    if(args()<5 | engine=="") engine=st_global("rdistnn_jackknife_engine")
    if(engine=="") engine="auto"
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    if(n<2 | m<1 | rows(signal)!=n | q<1 | hasmissing(signal) |
        hasmissing(ranks) | min(ranks)<1 | max(ranks)>n |
        any(ranks:!=floor(ranks)) | rows(table)!=n-1 |
        cols(table)!=1 | hasmissing(table) | (center!=0 & center!=1)) {
        errprintf("Rank-table deletions require finite aligned signals, rank permutations and n-1 weights.\n")
        exit(error(198))
    }
    for(j=1;j<=m;j++) {
        seen=J(n,1,0)
        for(i=1;i<=n;i++) {
            r=ranks[i,j]
            if(seen[r]) {
                errprintf("Each ranks column must be a permutation of 1..n.\n")
                exit(error(198))
            }
            seen[r]=1
        }
    }
    if(!any(engine:==("auto","plugin","mata","naive"))) {
        errprintf("rdistnn_jackknife_engine must be auto, plugin, mata or naive.\n")
        exit(error(198))
    }
    if(engine!="naive" & all(signal:==0)) return(J(n,m*q,0))
    if(center) {
        active=((table:!=0)\0):|(0\(table:!=0))
        loo=J(n,m*q,.)
        for(j=1;j<=m;j++) {
            z=rdistnn_center_signal(signal,active[ranks[,j]])
            loo[,((j-1)*q+1)..(j*q)]=rdistnn_loo_from_table(ranks[,j],z,table,0,engine)
        }
        return(loo)
    }
    if(engine=="naive") return(rdistnn_loo_naive_table(ranks,signal,table))
    plugin_mode=st_global("rdistnn_use_plugin")
    if(engine=="plugin" | (engine=="auto" & plugin_mode!="0")) {
        loo=rdistnn_loo_plugin(ranks,signal,table)
        if(rows(loo)==n & cols(loo)==m*q & !hasmissing(loo)) return(loo)
        if(engine=="plugin" | plugin_mode=="1") {
            errprintf("Compiled rank_loo unavailable or failed; install an updated plugin or allow Mata fallback.\n")
            exit(error(498))
        }
    }
    prefix=findexternal("rdistnn_loo_prefix()")
    if(prefix!=NULL) {
        loo=(*prefix)(ranks,signal,table)
        if(rows(loo)==n & cols(loo)==m*q & !hasmissing(loo)) return(loo)
    }
    if(engine=="mata") {
        errprintf("Mata prefix/suffix engine unavailable or produced nonfinite deletions.\n")
        exit(error(498))
    }
    return(rdistnn_loo_naive_table(ranks,signal,table))
}

// DNN adapter constructs a deletion table; the dispatcher owns computation.
real matrix rdistnn_loo_from_ranks(real matrix ranks, real matrix signal,
    real scalar k, real scalar s, | string scalar method, real scalar center,
    real colvector precomputed_table)
{
    real scalar n
    real colvector table
    if(args()<5 | method=="") method="exact"
    if(args()<6) center=0
    n=rows(ranks)
    if((method!="exact" & method!="asymp") | n<2 | missing(k) |
        missing(s) | k<1 | k!=floor(k) | s<k | s>=n | s!=floor(s)) {
        errprintf("rdistnn_loo_from_ranks: require 1 <= k <= s < n and exact/asymp weights.\n")
        exit(error(198))
    }
    if(args()<7) table=rdistnn_weights_from_ranks((1::(n-1)),k,s,method)
    else table=precomputed_table
    return(rdistnn_loo_from_table(ranks,signal,table,center))
}

// A common signal level has zero jackknife covariance. Remove it before
// accumulating deletions, and never restore it on covariance-only paths.
// Use only the union of deletion support to choose the reference, and zero
// unused rows before subtraction. A deleted earlier rank can bring the next
// original rank into support. Signed tables require nonzero support.
real matrix rdistnn_center_signal(real matrix signal, | real colvector active,
    real rowvector level)
{
    real scalar j
    real rowvector lo, hi
    real matrix z, used
    if(args()<2) active=J(rows(signal),1,1)
    z=J(rows(signal),cols(signal),0)
    level=J(1,cols(signal),0)
    if(!any(active)) return(z)
    used=select(signal,active)
    lo=colmin(used); hi=colmax(used)
    for(j=1;j<=cols(signal);j++) {
        if(lo[j]>0) level[j]=lo[j]
        else if(hi[j]<0) level[j]=hi[j]
    }
    z[selectindex(active),]=used:-level
    return(z)
}

end
