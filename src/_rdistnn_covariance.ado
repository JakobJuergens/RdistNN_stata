version 17
// Shared numeric covariance policy; estimator adapters own statistical factors.
capture mata: mata drop rdistnn_cov_symmetrize()
capture mata: mata drop rdistnn_covariance()
capture mata: mata drop rdistnn_cov_crossprod()
capture mata: mata drop rdistnn_cov_recover()

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

// Revisit sensitive exceptional entries at original scale. Keep the existing
// normalized result whenever direct products need overflow/underflow protection.
real matrix rdistnn_cov_recover(real matrix V, real matrix scaled,
    real matrix original, real matrix original_low, real rowvector exponent,
    real scalar factor, real scalar center)
{
    real scalar n,q,j,k,h,c,average,recovered
    real rowvector norm,prepared,bad,pair
    real colvector z,difference,back,a,b,value
    real matrix high,low
    if(!rows(original)) return(V)
    n=rows(original); q=cols(original)
    high=original
    low=rows(original_low) ? original_low : J(n,q,0)
    norm=sqrt(diagonal(scaled))'
    prepared=bad=J(1,q,0)
    for(j=1;j<q;j++) for(k=j+1;k<=q;k++) {
        if((exponent[j]==0 & exponent[k]==0) | V[j,j]>=. | V[k,k]>=.) continue
        if(abs(scaled[j,k])>sqrt(2.2204460492503131e-16)*norm[j]*norm[k]) continue
        pair=(j,k)
        for(h=1;h<=2;h++) {
            c=pair[h]
            if(prepared[c]) continue
            prepared[c]=1
            if(center) {
                z=original[,c]:-rdistnn_signal_level(original[,c])
                average=rdistnn_compensated_sum(z)/n
                if(average>=.) {
                    bad[c]=1
                    continue
                }
                difference=z:-average
                back=difference:-z
                low[,c]=(z:-(difference:-back)):-(average:+back)
                high[,c]=difference
            }
        }
        if(bad[j] | bad[k]) continue
        a=(high[,j] \ high[,j] \ low[,j] \ low[,j])
        b=(high[,k] \ low[,k] \ high[,k] \ low[,k])
        value=a:*b
        if(hasmissing(value) | any((a:!=0):&(b:!=0):&(abs(value):<2.2250738585072014e-308))) continue
        recovered=rdistnn_compensated_sum(value)
        if(recovered<. & factor*recovered<.) V[j,k]=V[k,j]=factor*recovered
    }
    return(V)
}

real matrix rdistnn_covariance(real matrix X, real scalar factor,
    real scalar center, | real scalar allow_failure)
{
    real scalar n,q,j,h,spread,extra,limit,power,half,value,original_center
    real rowvector exponent,invalid,varying,variance,means,amplitude,needs_scale
    real colvector difference,back
    real matrix Z,V,low,original,original_low,scaled
    if(args()<4) allow_failure=0
    n=rows(X); q=cols(X); Z=X
    exponent=J(1,q,0); invalid=colsum(Z:>=.):>0
    for(j=1;j<=q;j++) if(invalid[j]) Z[,j]=J(n,1,0)
    // The same conservative gate and arithmetic as R/utils-covariance.R.
    limit=sqrt(rdistnn_power2(1020)/n/max((1,factor)))
    amplitude=J(1,q,0)
    for(j=1;j<=q;j++) amplitude[j]=max(abs(Z[,j]))
    needs_scale=(amplitude:>=limit):|((amplitude:>0):&(amplitude:<sqrt(2.2250738585072014e-308)))
    original=original_low=J(0,0,0); original_center=center
    if(q>1 & any(needs_scale)) original=Z
    for(j=1;j<=q;j++) {
        if(needs_scale[j]) {
            exponent[j]=max((-1074,min((1022,floor(ln(amplitude[j])/ln(2))))))
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
            if(q>1 & !rows(original)) {
                original=Z; original_low=low; original_center=0
            }
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
        rdistnn_fail(498,"Covariance is outside the numeric range; rescale the signals or request point estimates only.")
    }
    scaled=rdistnn_cov_crossprod(Z,low)
    V=factor*scaled
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
    return(rdistnn_cov_recover(V,scaled,original,original_low,exponent,factor,original_center))
}
end
