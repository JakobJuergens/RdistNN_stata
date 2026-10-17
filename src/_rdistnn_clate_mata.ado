// Paired CLATE evaluator. The caller loads _rdistnn_mata first; reloading
// this extension drops only its own functions and preserves the base library.
capture mata: mata drop rdistnn_clate_fit()
capture mata: mata drop rdistnn_clate_diagnostics()

mata:

// Last five arguments are outputs. componentV is interleaved Y1,D1,Y2,D2.
// Predictions, ranks, metric and s remain fixed across paired deletions.
void rdistnn_clate_fit(real matrix ranks, real matrix signals,
    real scalar k, real scalar s, string scalar method,
    real scalar do_se, real scalar tol, real matrix components,
    real rowvector theta, real matrix V, real matrix componentV,
    real matrix status)
{
    real scalar n, m, i, j, h, d, den_scale, scale, variance, corr
    real colvector idx, kept, deleted_ranks, w, sd, good
    real matrix W, loo, centered, g, normalized, block

    n=rows(ranks); m=cols(ranks)
    if (n<1 | m<1 | rows(signals)!=n | cols(signals)!=2 |
        hasmissing(ranks) | hasmissing(signals) | tol<0 | tol>=. |
        k<1 | k!=floor(k) | s<k | s!=floor(s) | s>n |
        (do_se & s>=n)) {
        errprintf("rdistnn_clate_fit: invalid ranks, signals, tuning, or tolerance.\n")
        exit(error(198))
    }
    W=rdistnn_weights_from_rmat(ranks,k,s,method)
    components=W'*signals
    theta=J(1,m,.)
    V=J(m,m,.)
    componentV=J(2*m,2*m,.)
    status=J(m,5,0)
    status[,4]=J(m,1,1)
    den_scale=max(abs(signals[,2]))
    for(j=1;j<=m;j++) {
        d=components[j,2]
        status[j,2]=(d<0)
        if(d!=0 & d<. & den_scale>0 & tol>0) {
            status[j,3]=(abs(d)/den_scale<=tol)
        }
        if(hasmissing(components[j,])) status[j,1]=2
        else if(d==0) status[j,1]=1
        else {
            theta[j]=components[j,1]/d
            if(theta[j]>=.) status[j,1]=2
        }
    }
    if(!do_se) return

    idx=(1::n)
    loo=J(n,2*m,.)
    for(i=1;i<=n;i++) {
        kept=(idx:!=i)
        for(j=1;j<=m;j++) {
            deleted_ranks=select(ranks[,j],kept)
            deleted_ranks=deleted_ranks:-(deleted_ranks:>ranks[i,j])
            w=rdistnn_weights_from_ranks(deleted_ranks,k,s,method)
            loo[i,(2*j-1)..(2*j)]=cross(w,select(signals,kept))
        }
    }
    centered=loo:-mean(loo)
    componentV=(n-1)/n*cross(centered,centered)
    componentV=(componentV+componentV')/2
    normalized=J(n,m,.)
    sd=J(m,1,.)
    good=J(m,1,0)
    status[,4]=J(m,1,0)
    for(j=1;j<=m;j++) {
        block=centered[,(2*j-1)..(2*j)]
        if(hasmissing(block) |
            hasmissing(componentV[(2*j-1)..(2*j),(2*j-1)..(2*j)]) |
            (any(block[,1]:!=0) & componentV[2*j-1,2*j-1]==0) |
            (any(block[,2]:!=0) & componentV[2*j,2*j]==0)) {
            componentV[(2*j-1)..(2*j),]=J(2,2*m,.)
            componentV[,(2*j-1)..(2*j)]=J(2*m,2,.)
            status[j,4]=3
            continue
        }
        if(theta[j]>=.) continue
        // Linearize centered component deletions at the full-sample ratio.
        // Divide before combining to avoid squaring a tiny denominator.
        g=block[,1]/components[j,2]-theta[j]*(block[,2]/components[j,2])
        scale=max(abs(g))
        if(hasmissing(g) | scale>=.) {
            status[j,4]=3
            continue
        }
        if(scale==0) {
            normalized[,j]=J(n,1,0)
            sd[j]=0
            status[j,4]=2
        }
        else {
            normalized[,j]=g/scale
            sd[j]=scale*sqrt((n-1)/n*quadsum(normalized[,j]:^2))
        }
        variance=sd[j]^2
        if(variance>=. | (scale>0 & variance==0)) {
            status[j,4]=3
            continue
        }
        V[j,j]=variance
        good[j]=1
        status[j,5]=(variance>0)
    }
    for(j=1;j<=m;j++) if(good[j]) {
        for(h=j+1;h<=m;h++) if(good[h]) {
            if(sd[j]==0 | sd[h]==0) V[j,h]=0
            else {
                corr=cross(normalized[,j],normalized[,h]) /
                    sqrt(quadsum(normalized[,j]:^2)*quadsum(normalized[,h]:^2))
                V[j,h]=sd[j]*(corr*sd[h])
            }
            V[h,j]=V[j,h]
        }
    }
}

// Preserve CATE localization columns. Replace overlap by instrument arms
// and score concentration by |w*(phi_Y-theta*phi_D)|. Empty A/e mark the
// supplied-signal route, where overlap information is unavailable.
real matrix rdistnn_clate_diagnostics(real matrix diag, real matrix W,
    real matrix signals, real rowvector theta, real matrix A, real matrix e)
{
    real scalar n,m,j,total,mass1,mass0
    real matrix out
    real colvector w,armw,contrib,ord
    n=rows(W); m=cols(W)
    if(rows(diag)!=m | cols(diag)!=25 | rows(signals)!=n |
        cols(signals)!=2 | cols(theta)!=m |
        (rows(A)!=0 & (rows(A)!=n | cols(A)!=1)) |
        (rows(e)!=0 & (rows(e)!=n | cols(e)!=1))) {
        errprintf("rdistnn_clate_diagnostics: incompatible input dimensions.\n")
        exit(error(198))
    }
    out=diag
    out[,12..19]=J(m,8,.)
    out[,21..25]=J(m,5,.)
    out[,(12,16,21,24)]=J(m,4,0)
    for(j=1;j<=m;j++) {
        w=W[,j]
        if(rows(A)>0) {
            mass1=cross(A,w); mass0=cross(1:-A,w)
            out[j,12]=out[j,21]=1
            out[j,13..15]=(mass1,mass0,mass1/(mass1+mass0))
            armw=select(w,A:==1)
            if(rows(armw)>0 & quadsum(armw:^2)>0) {
                out[j,22]=quadsum(armw)^2/quadsum(armw:^2)
            }
            armw=select(w,A:==0)
            if(rows(armw)>0 & quadsum(armw:^2)>0) {
                out[j,23]=quadsum(armw)^2/quadsum(armw:^2)
            }
        }
        if(rows(e)>0) {
            out[j,16]=1
            out[j,17]=cross(e,w)/quadsum(w)
            out[j,18]=rdistnn_weighted_quantile(e,w,.1)
            out[j,19]=rdistnn_weighted_quantile(e,w,.9)
        }
        if(theta[j]<.) {
            contrib=abs(w:*(signals[,1]-theta[j]*signals[,2]))
            if(!hasmissing(contrib)) {
                out[j,24]=1
                total=quadsum(contrib)
                if(total>0 & total<.) {
                    ord=order((contrib,(1::n)),(-1,2))
                    out[j,25]=quadsum(contrib[ord[1..min((5,n))]])/total
                }
            }
        }
    }
    return(out)
}
end
