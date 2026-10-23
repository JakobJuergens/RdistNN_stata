*! version 0.1.0  2026-10-21
version 17
// Exact equality blocks. Highlevel callers supply one metric-distance column
// per target. Generic block primitives can also group explicit equality features.
foreach fn in distance_groups group_ties distance_warning distance_weights distance_loo duplicate_rows block_ranges block_weights block_add block_table_mean block_loo_direct block_loo_prefix block_loo_plugin block_loo {
    capture mata: mata drop rdistnn_`fn'()
}
mata:
// Group IDs follow sorted metric distance, independently of supplied ranks.
// Equality uses represented
// metric distances exactly. No coordinate or tolerance criterion is applied.
real matrix rdistnn_distance_groups(real matrix ranks, real matrix DX, | real scalar ordered)
{
    real scalar n,m,j,r,g
    real colvector ord
    real matrix groups
    n=rows(ranks); m=cols(ranks)
    if(rows(DX)!=n | cols(DX)!=m | hasmissing(DX)) exit(error(198))
    if(args()<3) ordered=0
    groups=J(n,m,.)
    for(j=1;j<=m;j++) {
        if(ordered) {
            ord=J(n,1,.)
            ord[ranks[,j]]=(1::n)
        }
        else ord=order((DX[,j],(1::n)),(1,2))
        g=1
        groups[ord[1],j]=g
        for(r=2;r<=n;r++) {
            if(DX[ord[r],j]!=DX[ord[r-1],j]) g++
            groups[ord[r],j]=g
        }
    }
    return(groups)
}

// Group IDs are contiguous 1..G, supplied by preparation or validated e(model_groups).
real rowvector rdistnn_group_ties(real matrix groups)
{
    return(colmax(groups):<rows(groups))
}

void rdistnn_distance_warning(real matrix groups)
{
    if(any(rdistnn_group_ties(groups))) {
        displayas("text")
        display("warning: Exact distance ties detected; weights are averaged within each target's equal-distance blocks, including jackknife deletions. The theory assumes a continuous distance distribution and does not establish inference for persistent ties.")
    }
}

real matrix rdistnn_distance_weights(real matrix ranks, real matrix W, real matrix DX, | real matrix prepared_groups)
{
    real scalar j
    real matrix out
    if(rows(DX)!=rows(ranks) | cols(DX)!=cols(ranks) | hasmissing(DX)) exit(error(198))
    out=W
    if(args()<4) prepared_groups=J(0,0,.)
    for(j=1;j<=cols(ranks);j++) out[,j]=rdistnn_block_weights(ranks[,j],W[,j],DX[,j],rows(prepared_groups) ? prepared_groups[,j] : J(0,1,.))
    return(out)
}

// The native ABI accepts one group column. Dispatch targets independently so
// each target retains its own exact distance groups and deletion support.
real matrix rdistnn_distance_loo(real matrix ranks, real matrix signal,
    real colvector table, real scalar center, string scalar engine, real matrix DX, | real matrix prepared_groups)
{
    real scalar j,q
    real matrix out
    if(rows(DX)!=rows(ranks) | cols(DX)!=cols(ranks) | hasmissing(DX)) exit(error(198))
    if(args()<7) prepared_groups=J(0,0,.)
    q=cols(signal); out=J(rows(ranks),cols(ranks)*q,.)
    for(j=1;j<=cols(ranks);j++) {
        if(rows(prepared_groups) ? max(prepared_groups[,j])<rows(ranks) : rdistnn_duplicate_rows(DX[,j])) {
            out[,((j-1)*q+1)..(j*q)]=rdistnn_block_loo(ranks[,j],signal,table,center,engine,DX[,j],rows(prepared_groups) ? prepared_groups[,j] : J(0,1,.))
        }
        else out[,((j-1)*q+1)..(j*q)]=rdistnn_loo_from_table(ranks[,j],signal,table,center,engine)
    }
    return(out)
}

real scalar rdistnn_duplicate_rows(real matrix X)
{
    real matrix sorted
    real scalar i
    if(rows(X)<2) return(0)
    sorted=X[order(X,(1..cols(X))),]
    for(i=2;i<=rows(X);i++) if(all(sorted[i,]:==sorted[i-1,])) return(1)
    return(0)
}

// Rows of the result are contiguous [first,last] rank intervals. Equality is
// checked on every supplied equality-feature column, never on signal values.
real matrix rdistnn_block_ranges(real colvector ranks, real matrix X)
{
    real colvector ord
    real scalar n,a,b,g
    real matrix out
    n=rows(ranks)
    if(rows(X)!=n | cols(X)<1 | hasmissing(X)) exit(error(198))
    ord=J(n,1,.); ord[ranks]=(1::n)
    out=J(n,2,.); a=1; g=0
    while(a<=n) {
        b=a
        while(b<n) {
            if(!all(X[ord[b+1],]:==X[ord[a],])) break
            b++
        }
        out[++g,]=(a,b)
        a=b+1
    }
    return(out[1..g,])
}

real matrix rdistnn_block_weights(real matrix ranks, real matrix W, real matrix X, | real colvector prepared_groups)
{
    real scalar j,g,a,b,mass
    real colvector ord,ix,rawranks
    real matrix blocks,out
    out=W
    if(rows(X)!=rows(ranks) | cols(X)<1 | hasmissing(X)) exit(error(198))
    if(args()<4) prepared_groups=J(0,1,.)
    if(rows(prepared_groups)) {
        if(max(prepared_groups)==rows(ranks)) return(out)
        ord=J(rows(ranks),1,.); ord[ranks[,1]]=(1::rows(ranks))
    }
    else {
        if(!rdistnn_duplicate_rows(X)) return(out)
        ord=order((X,(1::rows(X))),(1..(cols(X)+1)))
    }
    rawranks=J(rows(X),1,.); rawranks[ord]=(1::rows(X))
    blocks=rdistnn_block_ranges(rawranks,rows(prepared_groups) ? prepared_groups : X)
    for(j=1;j<=cols(ranks);j++) {
        for(g=1;g<=rows(blocks);g++) {
            a=blocks[g,1]; b=blocks[g,2]
            if(b==a) continue
            ix=ord[a..b]
            mass=rdistnn_sum(W[ix,j])/(b-a+1)
            // Sum first retains subnormal block mass; divide first only to
            // recover an overflowing signed/exceptional intermediate sum.
            if(missing(mass)) mass=rdistnn_sum(W[ix,j]:/(b-a+1))
            out[ix,j]=J(b-a+1,1,mass)
        }
    }
    return(out)
}

void rdistnn_block_add(real scalar high, real scalar low, real scalar value)
{
    real scalar next_high
    next_high=high+value
    if(next_high<.) {
        if(abs(high)>=abs(value)) low=low+((high-next_high)+value)
        else low=low+((value-next_high)+high)
    }
    else low=0
    high=next_high
}

// Inclusive 1-based table bounds; padded endpoints are used only by blocks
// that cannot occur before/after a particular deletion. Empty survivors give 0.
real scalar rdistnn_block_table_mean(real colvector table, real scalar first, real scalar last)
{
    real colvector values
    real scalar out,c,a,b
    if(last<first) return(0)
    c=last-first+1
    a=max((1,first)); b=min((rows(table),last))
    values=a<=b ? table[a..b] : J(0,1,0)
    out=rdistnn_sum(values)/c
    if(missing(out)) out=rdistnn_sum(values:/c)
    return(out)
}

// Explicit regrouping is also the bounded fallback for supplied ranks whose
// equality groups interleave or reorder members. This path takes O(n^2*q).
real matrix rdistnn_block_loo_direct(real matrix z, real colvector table, real colvector groups,
    | real colvector observation_ids)
{
    real scalar n,q,i,r,g,high,low,value,rank
    real colvector fullcount,count,leading,correction,weights,kept,bad,observation_order,mass_order
    real matrix out
    n=rows(z); q=cols(z); out=J(n,q,0)
    if(args()<4) observation_ids=(1::n)
    observation_order=order(observation_ids,1)
    fullcount=J(n,1,0)
    for(r=1;r<=n;r++) fullcount[groups[r]]=fullcount[groups[r]]+1
    for(i=1;i<=n;i++) {
        kept=selectindex((1::n):!=i)
        mass_order=select(observation_order,observation_order:!=i)
        count=fullcount; count[groups[i]]=count[groups[i]]-1
        leading=correction=J(n,1,0)
        // Public block weights reduce assigned entries in original row order.
        // Keep that order after the outer dispatcher has sorted z by rank.
        for(r=1;r<=n-1;r++) {
            rank=mass_order[r]; g=groups[rank]; high=leading[g]; low=correction[g]
            rdistnn_block_add(high,low,table[rank-(rank>i)])
            leading[g]=high; correction[g]=low
        }
        weights=J(n,1,0)
        for(g=1;g<=n;g++) if(count[g]) weights[g]=(leading[g]+correction[g])/count[g]
        bad=(weights:>=.)
        if(any(bad)) {
            leading=correction=J(n,1,0)
            for(r=1;r<=n-1;r++) {
                rank=mass_order[r]; g=groups[rank]
                if(!bad[g]) continue
                high=leading[g]; low=correction[g]; value=table[rank-(rank>i)]/count[g]
                rdistnn_block_add(high,low,value)
                leading[g]=high; correction[g]=low
            }
            weights[selectindex(bad)]=leading[selectindex(bad)]+correction[selectindex(bad)]
        }
        out[i,]=rdistnn_weighted_signals(z[kept,],weights[groups[kept]])'
    }
    return(out)
}

// Original-signal prefix/postfix states. Other blocks use their actual
// unshifted/shifted averaged weights; the deleted block uses survivor weights.
// The four high/low states are rounded only after their final common merge.
real matrix rdistnn_block_loo_prefix(real matrix z, real colvector table, real colvector groups,
    | real colvector observation_ids)
{
    real scalar n,q,a,b,g,G,r,h,high,low,right_high,right_low,value
    real colvector starts,ends,seen,before,after,survivor,terms
    real matrix ph,pl,sh,sl,lh,ll,out
    n=rows(z); q=cols(z); starts=ends=J(n,1,.); seen=J(n,1,0)
    if(args()<4) observation_ids=(1::n)
    G=0; a=1
    while(a<=n) {
        g=groups[a]
        if(seen[g]) return(rdistnn_block_loo_direct(z,table,groups,observation_ids))
        seen[g]=1; b=a
        while(b<n) {
            if(groups[b+1]!=g) break
            if(observation_ids[b+1]<observation_ids[b]) {
                return(rdistnn_block_loo_direct(z,table,groups,observation_ids))
            }
            b++
        }
        starts[++G]=a; ends[G]=b; a=b+1
    }
    starts=starts[1..G]; ends=ends[1..G]
    before=after=survivor=J(G,1,0)
    for(g=1;g<=G;g++) {
        before[g]=rdistnn_block_table_mean(table,starts[g],ends[g])
        after[g]=rdistnn_block_table_mean(table,starts[g]-1,ends[g]-1)
        survivor[g]=rdistnn_block_table_mean(table,starts[g],ends[g]-1)
    }
    ph=pl=sh=sl=J(G+1,q,0); lh=ll=J(n,q,0); out=J(n,q,0)
    for(h=1;h<=q;h++) {
        high=low=0
        for(g=1;g<=G;g++) {
            ph[g,h]=high; pl[g,h]=low
            for(r=starts[g];r<=ends[g];r++) {
                value=z[r,h]*before[g]; rdistnn_block_add(high,low,value)
            }
        }
        high=low=0
        for(g=G;g>=1;g--) {
            sh[g+1,h]=high; sl[g+1,h]=low
            for(r=ends[g];r>=starts[g];r--) {
                value=z[r,h]*after[g]; rdistnn_block_add(high,low,value)
            }
        }
        for(g=1;g<=G;g++) {
            high=low=0
            for(r=starts[g];r<=ends[g];r++) {
                lh[r,h]=high; ll[r,h]=low
                value=z[r,h]*survivor[g]; rdistnn_block_add(high,low,value)
            }
            right_high=right_low=0
            for(r=ends[g];r>=starts[g];r--) {
                terms=(ph[g,h]\sh[g+1,h]\lh[r,h]\right_high\pl[g,h]\sl[g+1,h]\ll[r,h]\right_low)
                out[r,h]=rdistnn_compensated_sum(terms)
                value=z[r,h]*survivor[g]; rdistnn_block_add(right_high,right_low,value)
            }
        }
    }
    return(out)
}

real matrix rdistnn_block_loo_plugin(real matrix ranks, real matrix z, real colvector table, real colvector groups)
{
    if(st_global("rdistnn_plugin_ok")!="1" | !rdistnn_native_packed_api()) return(J(0,0,.))
    return(rdistnn_native_loo(ranks,z,table,groups))
}

// Center on the union of deletion support, expanded to every member of each
// intersecting equality group. Each selected engine computes original-signal
// deletions itself; forcing native cannot run a dummy rank_loo calculation.
real matrix rdistnn_block_loo(real matrix ranks, real matrix signal,
    real colvector table, real scalar center, string scalar engine, real matrix X, | real colvector prepared_groups)
{
    real scalar n,q,j,g,a,b,h
    real colvector ord,active,raword,rawranks,positions,groups
    real rowvector level
    real matrix out,blocks,z,values
    string scalar plugin_mode
    n=rows(ranks); q=cols(signal)
    if(rows(X)!=n | cols(X)<1 | hasmissing(X)) exit(error(198))
    if(args()<7) prepared_groups=J(0,1,.)
    if(rows(prepared_groups)) {
        raword=J(n,1,.); raword[ranks[,1]]=(1::n)
        blocks=rdistnn_block_ranges(ranks[,1],prepared_groups)
        groups=prepared_groups
    }
    else {
        raword=order((X,(1::rows(X))),(1..(cols(X)+1)))
        rawranks=J(n,1,.); rawranks[raword]=(1::n)
        blocks=rdistnn_block_ranges(rawranks,X); groups=J(n,1,.)
        for(g=1;g<=rows(blocks);g++) groups[raword[blocks[g,1]..blocks[g,2]]]=J(blocks[g,2]-blocks[g,1]+1,1,g)
    }
    plugin_mode=st_global("rdistnn_use_plugin")
    out=J(n,cols(ranks)*q,.)
    for(j=1;j<=cols(ranks);j++) {
        ord=J(n,1,.); ord[ranks[,j]]=(1::n)
        active=((table:!=0)\0):|(0\(table:!=0))
        for(g=1;g<=rows(blocks);g++) {
            a=blocks[g,1]; b=blocks[g,2]
            positions=ranks[raword[a..b],j]
            if(any(active[positions])) active[positions]=J(b-a+1,1,1)
        }
        z=rdistnn_center_signal(signal[ord,],active,level,table)
        values=J(0,0,.)
        if(engine=="plugin" | (engine=="auto" & plugin_mode!="0")) {
            // Retain original observation identities through native transport;
            // noncanonical blocks use them to reproduce public block masses.
            values=rdistnn_block_loo_plugin(ranks[,j],z[ranks[,j],],table,groups)
            if(rows(values)==n & cols(values)==q & !hasmissing(values)) {
                values=values[ord,]
                st_global("rdistnn_last_jackknife_engine","plugin")
            }
            else {
                values=J(0,0,.)
                if(engine=="plugin" | plugin_mode=="1") {
                    errprintf("Compiled block deletion API is unavailable; install the current plugin or set rdistnn_jackknife_engine mata.\n")
                    exit(error(498))
                }
            }
        }
        if(!rows(values)) {
            if(engine=="naive") values=rdistnn_block_loo_direct(z,table,groups[ord],ord)
            else values=rdistnn_block_loo_prefix(z,table,groups[ord],ord)
            if(hasmissing(values)) {
                errprintf("Original-signal block deletions exceed the numeric range.\n")
                exit(error(498))
            }
            st_global("rdistnn_last_jackknife_engine",engine=="naive" ? "naive" : "mata")
        }
        if(!center) values=values:+level*rdistnn_sum(table)
        out[ord,((j-1)*q+1)..(j*q)]=values
    }
    return(out)
}
end
