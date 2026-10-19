version 17
// Shared numeric covariance policy; estimator adapters own statistical factors.
capture mata: mata drop rdistnn_cov_symmetrize()
capture mata: mata drop rdistnn_covariance()
capture mata: mata drop rdistnn_cov_crossprod()

mata:

real matrix rdistnn_cov_symmetrize(real matrix V)
{
    real scalar i,j,a,b,mid
    for(i=1;i<rows(V);i++) for(j=i+1;j<=cols(V);j++) {
        a=V[i,j]; b=V[j,i]
        if(a==b) mid=a
        else if((a>=0)==(b>=0)) mid=a+(b-a)/2
        else mid=a/2+b/2
        V[i,j]=V[j,i]=mid
    }
    return(V)
}

real matrix rdistnn_cov_crossprod(real matrix Z, real matrix low)
{
    real matrix V
    real rowvector norm
    real colvector value
    real scalar j,k
    V=cross(Z,Z)
    norm=sqrt(diagonal(V))'
    for(j=1;j<cols(Z);j++) for(k=j+1;k<=cols(Z);k++) {
        if(abs(V[j,k])<=sqrt(2.2204460492503131e-16)*norm[j]*norm[k]) {
            value=Z[,j]:*Z[,k]
            if(rows(low)) value=(value \ Z[,j]:*low[,k] \ low[,j]:*Z[,k] \ low[,j]:*low[,k])
            V[j,k]=V[k,j]=rdistnn_compensated_sum(value)
        }
    }
    return(V)
}

real matrix rdistnn_covariance(real matrix X, real scalar factor,
    real scalar center, | real scalar allow_failure)
{
    real scalar n,q,j,h,amplitude,spread,extra,limit,power,half,value
    real rowvector exponent,invalid,varying,variance,means
    real colvector difference,back
    real matrix Z,V,low
    if(args()<4) allow_failure=0
    n=rows(X); q=cols(X); Z=X
    exponent=J(1,q,0); invalid=colsum(Z:>=.):>0
    for(j=1;j<=q;j++) if(invalid[j]) Z[,j]=J(n,1,0)
    // The same conservative gate and arithmetic as R/utils-covariance.R.
    limit=sqrt(rdistnn_power2(1020)/n/max((1,factor)))
    for(j=1;j<=q;j++) {
        amplitude=max(abs(Z[,j]))
        if(amplitude>=limit | (amplitude>0 & amplitude<sqrt(2.2250738585072014e-308))) {
            exponent[j]=max((-1074,min((1022,floor(ln(amplitude)/ln(2))))))
            Z[,j]=Z[,j]/rdistnn_power2(exponent[j])
        }
        if(center) Z[,j]=Z[,j]:-rdistnn_signal_level(Z[,j])
    }
    low=J(0,0,0)
    if(center) {
        means=mean(Z)
        if(q>1) low=J(n,q,0)
        for(j=1;j<=q;j++) {
            difference=Z[,j]:-means[j]
            if(q>1) {
                back=difference:-Z[,j]
                low[,j]=(Z[,j]:-(difference:-back)):-(means[j]:+back)
            }
            Z[,j]=difference
        }
    }
    varying=colsum(Z:!=0):>0
    for(j=1;j<=q;j++) {
        spread=max(abs(Z[,j]))
        if(spread>0 & spread<sqrt(2.2250738585072014e-308)) {
            extra=max((-1074,min((1022,floor(ln(spread)/ln(2))))))
            Z[,j]=Z[,j]/rdistnn_power2(extra)
            if(rows(low)) low[,j]=low[,j]/rdistnn_power2(extra)
            exponent[j]=exponent[j]+extra
        }
    }
    variance=factor*quadcolsum(Z:*Z)
    for(j=1;j<=q;j++) {
        power=2*exponent[j]; half=trunc(power/2)
        variance[j]=(variance[j]*rdistnn_power2(half))*rdistnn_power2(power-half)
    }
    invalid=invalid :| (variance:>=.) :| (variance:<0) :| (varying :& (variance:==0))
    if(any(invalid) & !allow_failure) {
        errprintf("Covariance is outside the numeric range; rescale the signals or request point estimates only.\n")
        exit(error(498))
    }
    V=factor*rdistnn_cov_crossprod(Z,low)
    for(j=1;j<=q;j++) for(h=j;h<=q;h++) {
        power=exponent[j]+exponent[h]; half=trunc(power/2)
        value=(V[j,h]*rdistnn_power2(half))*rdistnn_power2(power-half)
        V[j,h]=V[h,j]=value
    }
    for(j=1;j<=q;j++) {
        V[j,j]=variance[j]
        if(invalid[j]) {
            V[j,]=J(1,q,.)
            V[,j]=J(q,1,.)
        }
    }
    return(V)
}
end
