version 17
// Paired CLATE evaluator. The caller loads _rdistnn_mata first; reloading
// this extension drops only its own functions and preserves the base library.
capture mata: mata drop rdistnn_clate_fit()

mata:

// Last five arguments are outputs. componentV is interleaved Y1,D1,Y2,D2.
// Predictions, ranks, metric and s remain fixed across paired deletions.
void rdistnn_clate_fit(real matrix ranks, real matrix signals,
    real scalar k, real scalar s, string scalar method,
    real scalar do_se, real scalar tol, real matrix components,
    real rowvector theta, real matrix V, real matrix componentV,
    real matrix status, | real matrix DX, real matrix prepared_groups)
{
    real scalar n, m, j, h, d, den_scale, scale, variance, corr, omega_quarter
    real colvector sd, good
    real matrix W, loo, centered, g, normalized, block

    if(args()<14) prepared_groups=J(0,0,.)
    n=rows(ranks); m=cols(ranks)
    if (n<1 | m<1 | rows(signals)!=n | cols(signals)!=2 |
        hasmissing(ranks) | hasmissing(signals) | tol<0 | tol>=. |
        k<1 | k!=floor(k) | s<k | s!=floor(s) | s>n |
        (do_se & s>=n)) {
        errprintf("rdistnn_clate_fit: invalid ranks, signals, tuning, or tolerance.\n")
        exit(error(198))
    }
    if(args()>=13 & rows(DX)>0) W=rdistnn_weights_from_rmat(ranks,k,s,method,DX,prepared_groups)
    else W=rdistnn_weights_from_rmat(ranks,k,s,method)
    components=rdistnn_weighted_signals(signals,W)'
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

    if(args()>=13 & rows(DX)>0) {
        loo=rdistnn_loo_from_ranks(ranks,signals,k,s,method,1,
            rdistnn_weights_from_ranks((1::(n-1)),k,s,method),DX,prepared_groups)
    }
    else loo=rdistnn_loo_from_ranks(ranks,signals,k,s,method,1)
    centered=loo:-mean(loo)
    componentV=rdistnn_covariance(loo,(n-1)/n,1,1)
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
            sd[j]=sqrt(rdistnn_covariance(g,(n-1)/n,0,1)[1,1])
        }
        variance=sd[j]^2
        if(variance>=. | (scale>0 & variance==0)) {
            status[j,4]=3
            continue
        }
        // R also requires finite neff*variance. Stata's nonmissing range is
        // smaller than R's IEEE range, so compare one quarter of that value
        // with DBL_MAX/4 instead of forming an unsupported omega directly.
        omega_quarter=variance*((1/quadsum(W[,j]:^2))/4)
        if(omega_quarter>4.4942328371557893e307) {
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


end
