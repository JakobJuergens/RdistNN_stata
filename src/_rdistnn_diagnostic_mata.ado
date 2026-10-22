version 17
// Matrix-input diagnostics shared by fresh commands and saved numeric models.
// Parsing, labels, display and estimator inference policy stay with callers.
capture mata: mata drop rdistnn_grid_weight_diagnostics()
capture mata: mata drop rdistnn_score_share()
capture mata: mata drop rdistnn_diagnostics()
capture mata: mata drop rdistnn_clate_diagnostics()
mata:
// Retain the five-column candidate summary, taking concentration extrema
// across actual target-specific weights after duplicate-block averaging.
real rowvector rdistnn_grid_weight_diagnostics(real matrix W)
{
    real scalar j,n,top
    real rowvector ess,mx,share
    real colvector w
    n=rows(W); top=min((5,n))
    ess=1:/colsum(W:^2); mx=colmax(abs(W)); share=J(1,cols(W),.)
    for(j=1;j<=cols(W);j++) {
        w=sort(abs(W[,j]),-1)
        share[j]=rdistnn_sum(w[1..top])/rdistnn_sum(w)
    }
    return((min(ess),max(mx),top,max(share),any(ess:>=.9*n)))
}

// Concentration of |w*(left-right)|, with right already rounded by its
// caller. Recover only finite-operand range failures; no overflowing product
// or uncertain cancellation magnitude is inferred.
real rowvector rdistnn_score_share(real colvector w, real colvector left,
    real colvector right, | real colvector contribution)
{
    real scalar n,i,total,top,power,value,remaining,step
    real colvector active,ids,weights,l,r,raw,contrib,mantissa,exponent,direction,ord
    real rowvector a,b,c,parts
    if(rows(w)!=rows(left) | rows(w)!=rows(right)) exit(error(198))
    if(args()>=4) contribution=J(rows(w),1,.)
    if(hasmissing(w) | any(w:<0)) return((0,.))
    active=selectindex(w:>0)
    if(args()>=4) contribution=J(rows(w),1,0)
    if(!length(active)) return((1,.))
    ids=active
    weights=w[active]; l=left[active]; r=right[active]; n=rows(weights)
    raw=weights:*(l-r)
    if(args()>=4) contribution[ids]=raw
    if(hasmissing(l) | hasmissing(r)) return((0,.))
    contrib=abs(raw); total=quadsum(contrib,1)
    if(!hasmissing(contrib) & total<.) {
        if(total==0) return((1,.))
        ord=order((contrib,(1::n)),(-1,2))
        return((1,quadsum(contrib[ord[1..min((5,n))]],1)/total))
    }
    mantissa=exponent=J(n,1,0); direction=J(n,1,1)
    for(i=1;i<=n;i++) {
        if(contrib[i]<.) {
            if(contrib[i]==0) continue
            parts=rdistnn_binary_parts(contrib[i])
            mantissa[i]=parts[1]; exponent[i]=parts[2]
        }
        else {
            // Overflowed subtraction of two finite operands is an opposite-
            // sign, noncancelling case. Keep tiny weights in exponent space.
            a=rdistnn_binary_parts(l[i]); b=rdistnn_binary_parts(r[i])
            power=max((a[2],b[2]))
            if(l[i]==0) power=b[2]
            if(r[i]==0) power=a[2]
            value=l[i]/rdistnn_power2(power)-r[i]/rdistnn_power2(power)
            c=rdistnn_binary_parts(weights[i])
            value=value*c[1]
            direction[i]=value<0 ? -1 : 1
            value=abs(value)
            if(missing(value) | value==0) return((0,.))
            parts=rdistnn_binary_parts(value)
            mantissa[i]=parts[1]; exponent[i]=power+c[2]+parts[2]
        }
    }
    if(args()>=4) {
        for(i=1;i<=n;i++) if(raw[i]>=.) {
            value=direction[i]*mantissa[i]; remaining=exponent[i]
            // Keep the first negative step normal. A -1074 first step could
            // round a subnormal twice before restoring the remaining power.
            while(remaining!=0 & value<.) {
                step=max((-1022,min((1022,remaining))))
                value=value*rdistnn_power2(step)
                remaining=remaining-step
            }
            contribution[ids[i]]=value
        }
    }
    active=selectindex(mantissa:>0)
    if(!length(active)) return((1,.))
    top=max(exponent[active]); contrib=J(n,1,0)
    for(i=1;i<=length(active);i++) {
        power=exponent[active[i]]-top
        if(power>=(-1074)) contrib[active[i]]=mantissa[active[i]]*rdistnn_power2(power)
    }
    total=quadsum(contrib,1)
    ord=order((contrib,(1::n)),(-1,2))
    return((1,quadsum(contrib[ord[1..min((5,n))]],1)/total))
}

real matrix rdistnn_diagnostics(string scalar cmd, real matrix X,
    real matrix X0, real matrix signals, real matrix W, real matrix DX,
    real scalar k, real scalar s, real rowvector b, real matrix A, real matrix E)
{
    real scalar n,m,j,total,ybar
    real matrix out,pair
    real colvector w,d,ord,sorted
    n=rows(X); m=rows(X0)
    out=J(m,25,.)
    if(cmd=="np_dnn") out=J(m,17,.)
    if(cmd=="plqr_dnn") out=J(m,11,.)
    for(j=1;j<=m;j++) {
        w=W[,j]; d=DX[,j]; total=quadsum(w)
        ord=order((w,(1::n)),(-1,2)); sorted=sort(d,1)
        out[j,1..11]=(1/quadsum(w:^2),max(w),w[ord[1]],
            quadsum(w[ord[1..min((5,n))]]),quadsum(w[ord[1..min((10,n))]]),
            sum(abs(w):>sqrt(epsilon(1))),1,sorted[min((k,n))],
            sorted[min((s,n))],cross(d,w)/total,rdistnn_weighted_quantile(d,w,.5))
        if(cmd=="np_dnn") {
            ybar=rdistnn_weighted_signals(signals,w)/total
            out[j,12..17]=(1,ybar,rdistnn_weighted_sd(signals,w),
                rdistnn_weighted_quantile(signals,w,.1),
                rdistnn_weighted_quantile(signals,w,.5),
                rdistnn_weighted_quantile(signals,w,.9))
        }
        else if(cmd!="plqr_dnn") {
            out[j,20]=rdistnn_weighted_quantile(rdistnn_physical_distances(X,X0[j,]),w,.9)
        }
    }
    if(cmd=="cate_dnn" | cmd=="clate_dnn") {
        pair=signals
        if(cmd=="cate_dnn") pair=(signals,J(n,1,1))
        out=rdistnn_clate_diagnostics(out,W,pair,b,A,E)
    }
    return(out)
}

// Preserve CATE localization columns. Replace overlap by instrument arms
// and score concentration by |w*(phi_Y-theta*phi_D)|. Empty A/e mark the
// supplied-signal route, where overlap information is unavailable.
real matrix rdistnn_clate_diagnostics(real matrix diag, real matrix W,
    real matrix signals, real rowvector theta, real matrix A, real matrix e)
{
    real scalar n,m,j,total,mass1,mass0
    real matrix out
    real colvector w,armw,active
    real rowvector share
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
            active=selectindex(w:>0); share=(1,.)
            if(length(active)) share=rdistnn_score_share(w[active],
                signals[active,1],theta[j]*signals[active,2])
            out[j,24..25]=share
        }
    }
    return(out)
}
end
