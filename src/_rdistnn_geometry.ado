version 17
// Fitted Mahalanobis whitening from the centered training design, not its Gram.
capture mata: mata drop rdistnn_mahal_geometry()
capture mata: mata drop rdistnn_mahal_project()

mata:

// Packed geometry: scale/rank metadata, raw safe origin, then the
// p x p whitening factor (unused columns zero). The origin can be zero
// when an observed origin would lose represented differences.
real matrix rdistnn_mahal_geometry(real matrix X)
{
    real scalar n,p,amplitude,scale,tolerance,j,power,spread,global_power
    real rowvector origin,metadata,column_power,spread_power
    real colvector values,keep
    real matrix Z,U,Vt,V,B
    n=rows(X); p=cols(X)
    if(n<2) {
        errprintf("Mahalanobis distance requires at least two training observations.\n")
        exit(error(2001))
    }
    if(p<1) {
        errprintf("Mahalanobis distance requires at least one covariate.\n")
        exit(error(198))
    }
    // Retain an exact-subtraction-safe origin for projection and prediction;
    // keep the centered SVD arithmetic and numerical rank policy unchanged.
    origin=rdistnn_signal_level(X); Z=J(n,p,0)
    column_power=J(1,p,0); spread_power=J(1,p,.)
    // Choose a common scale from centered variation, not uncentered levels.
    // A large constant must not erase a small but resolvable varying column.
    for(j=1;j<=p;j++) {
        // A computed finite negative magnitude can exceed Stata's positive
        // finite ceiling: abs() is then missing and must not be skipped.
        if(hasmissing(abs(X[,j]))) power=1022
        else {
            amplitude=max(abs(X[,j]))
            power=amplitude==0 ? 0 : max((-1074,min((1022,floor(ln(amplitude)/ln(2))))))
        }
        Z[,j]=X[,j]/rdistnn_power2(power)
        Z[,j]=Z[,j]:-Z[1,j]; Z[,j]=Z[,j]:-mean(Z[,j])
        column_power[j]=power
        spread=max(abs(Z[,j]))
        if(spread>0) spread_power[j]=power+floor(ln(spread)/ln(2))
    }
    global_power=0
    if(any(spread_power:<.)) global_power=max((-1074,min((1022,max(select(spread_power,spread_power:<.))))))
    scale=rdistnn_power2(global_power)
    for(j=1;j<=p;j++) if(spread_power[j]<.) {
        power=column_power[j]-global_power
        // Values far below the common rank cutoff can round to zero.
        if(power < -1074) Z[,j]=J(n,1,0)
        else Z[,j]=Z[,j]*rdistnn_power2(power)
    }
    // Mata's economy SVD requires at least as many rows as columns.
    if(n>=p) {
        svd(Z,U,values,Vt)
        V=Vt'
    }
    else svd(Z',V,values,Vt)
    tolerance=max((n,p))*2.2204460492503131e-16*max(values)
    keep=selectindex(values:>tolerance)
    B=J(p,p,0)
    for(j=1;j<=length(keep);j++) B[,j]=V[,keep[j]]/(values[keep[j]]/sqrt(n-1))
    for(j=1;j<=p;j++) if(spread_power[j]>=.) B[j,]=J(1,p,0)
    metadata=J(1,p,0); metadata[1]=scale
    if(p>=2) metadata[2]=length(keep)
    return(metadata \ origin \ B)
}

real matrix rdistnn_mahal_project(real matrix X, real matrix geometry)
{
    real scalar j,c,origin,i,amplitude,power,remaining,step,finite_difference
    real colvector active,recover,difference,nonzero,exceptional,bad
    real rowvector value,row_difference,normalized
    real matrix Z,B
    B=geometry[|3,1 \ rows(geometry),cols(geometry)|]
    if(all(B:==0)) return(J(rows(X),cols(X),0))
    active=selectindex(rowsum(abs(B)):>0)
    Z=X[,active]
    exceptional=J(rows(X),1,0)
    for(j=1;j<=length(active);j++) {
        c=active[j]; origin=geometry[2,c]
        difference=Z[,j]:-origin
        nonzero=difference:!=0
        recover=selectindex((difference:>=.):&(geometry[1,1]>1))
        difference=difference/geometry[1,1]
        if(length(recover)) difference[recover]=Z[recover,j]/geometry[1,1]:-origin/geometry[1,1]
        exceptional=exceptional :| (difference:>=.) :|
            (nonzero :& (abs(difference):<2.2250738585072014e-308))
        Z[,j]=difference
    }
    bad=selectindex(exceptional)
    if(length(bad)) Z[bad,]=J(length(bad),cols(Z),0)
    // Reduce each projection independently of the matrix/vector BLAS shape.
    Z=rdistnn_weighted_signals(Z',B[active,])
    for(j=1;j<=length(bad);j++) {
        i=bad[j]
        row_difference=X[i,active]:-geometry[2,active]
        finite_difference=!hasmissing(row_difference)
        amplitude=finite_difference ? max(abs(row_difference)) : max(abs((X[i,active],geometry[2,active])))
        power=amplitude==0 ? 0 : max((-1074,min((1022,floor(ln(amplitude)/ln(2))))))
        if(finite_difference) normalized=row_difference/rdistnn_power2(power)
        else normalized=X[i,active]/rdistnn_power2(power):-geometry[2,active]/rdistnn_power2(power)
        value=rdistnn_weighted_signals(normalized',B[active,])
        remaining=power-floor(ln(geometry[1,1])/ln(2)+.5)
        while(remaining!=0) {
            step=max((-1074,min((1022,remaining))))
            value=value*rdistnn_power2(step)
            remaining=remaining-step
        }
        Z[i,]=value
    }
    if(hasmissing(Z)) {
        errprintf("Mahalanobis projection is outside the numeric range; rescale covariates and targets.\n")
        exit(error(498))
    }
    return(Z)
}
end
