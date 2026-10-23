*! version 0.1.0  2026-10-19
version 17
// Shared finite double-precision arithmetic; estimators own statistical policy.
foreach fn in signal_level sum compensated_sum weighted_cdf_index weighted_sd weighted_signals accumulation_add accumulation_merge accumulation_deletions {
    capture mata: mata drop rdistnn_`fn'()
}
capture mata: mata drop rdistnn_scale_sum()
capture mata: mata drop rdistnn_power2()
capture mata: mata drop rdistnn_binary_scale()
capture mata: mata drop rdistnn_binary_parts()
capture mata: mata drop rdistnn_weighted_quantile()
mata:
real rowvector rdistnn_signal_level(real matrix signal, | real matrix weights, real scalar weight_mass)
{
    real scalar j,amplitude,weight_amplitude
    real rowvector lo,hi,level
    real colvector difference,back,residual
    if(!rows(signal)) return(J(1,cols(signal),0))
    lo=colmin(signal); hi=colmax(signal)
    level=signal[1,]
    for(j=1;j<=cols(signal);j++) {
        if(abs(hi[j]-level[j])>=8.9884656743115795e307 |
           abs(lo[j]-level[j])>=8.9884656743115795e307) {
            level[j]=0
            if(lo[j]>0) level[j]=lo[j]
            else if(hi[j]<0) level[j]=hi[j]
        }
    }
    weight_amplitude=0
    if(args()>=2) {
        if(length(weights)) weight_amplitude=max(vec(abs(weights)))
        if(args()<3) weight_mass=quadsum(weights,1)
    }
    for(j=1;j<=cols(signal);j++) {
        // A subnormal origin can erase an unrepresentable weighted residual
        // before restoration. Match the direct separately-rounded products.
        amplitude=max((abs(lo[j]),abs(hi[j])))
        if(amplitude<2.2250738585072014e-308) level[j]=0
        if(args()>=2) {
            if(weight_amplitude*amplitude<2.2250738585072014e-308) level[j]=0
            if(weight_mass!=0) {
                if(abs(weight_mass)<2.2250738585072014e-308 |
                    abs(weight_mass)*amplitude<2.2250738585072014e-308) level[j]=0
            }
        }
    }
    for(j=1;j<=cols(signal);j++) {
        difference=signal[,j]:-level[j]
        back=difference:-signal[,j]
        residual=(signal[,j]:-(difference:-back)):-(level[j]:+back)
        if(any(residual:!=0)) level[j]=0
    }
    return(level)
}

real scalar rdistnn_sum(real matrix x)
{
    // Native quad reductions match the canonical compensated reduction on
    // shared fixtures. Explicitly propagate missing values, never omit them.
    return(quadsum(x,1))
}

// Canonical serial double reduction for decisions at a CDF boundary.
real scalar rdistnn_compensated_sum(real colvector x)
{
    real scalar i,high,low,next_high,value
    high=low=0
    for(i=1;i<=rows(x);i++) {
        value=x[i]; next_high=high+value
        if(next_high<.) {
            if(abs(high)>=abs(value)) low=low+((high-next_high)+value)
            else low=low+((value-next_high)+high)
        }
        else low=0
        high=next_high
    }
    return(high+low)
}

real scalar rdistnn_weighted_cdf_index(real colvector w, real scalar prob)
{
    real scalar mass,i,high,low,next_high,value
    mass=rdistnn_compensated_sum(w)
    if(missing(mass) | mass<=0) return(.)
    high=low=0
    for(i=1;i<=rows(w);i++) {
        value=w[i]; next_high=high+value
        if(abs(high)>=abs(value)) low=low+((high-next_high)+value)
        else low=low+((value-next_high)+high)
        high=next_high
        if((high+low)/mass>=prob) return(i)
    }
    return(rows(w))
}

// Population SD as a norm; the variance need not be representable. Match the
// R helper's binary scale, observed origin and separately rounded products.
real scalar rdistnn_weighted_sd(real colvector x, real colvector w)
{
    real scalar mass,amplitude,scale,center,spread,norm
    real colvector ix,z,weights,root_weights,d
    ix=selectindex(w:>0)
    if(any(w:<0) | !length(ix)) return(.)
    z=x[ix]; weights=w[ix]
    mass=rdistnn_compensated_sum(weights)
    if(missing(mass) | mass<=0) return(.)
    // Shift before normalization, with the common exact-subtraction guard.
    // A tiny-weight extreme cannot erase the other outcomes' local spread.
    z=z:-rdistnn_signal_level(z)
    if(hasmissing(abs(z))) {
        // A finite computed negative value can exceed Stata's positive
        // missing boundary in magnitude. Match R's capped binary scale.
        scale=rdistnn_power2(1022)
    }
    else {
        amplitude=max(abs(z))
        if(amplitude==0) return(0)
        scale=rdistnn_binary_scale(amplitude)
    }
    z=z:/scale
    // A normalized weight can underflow while its root and the SD remain
    // representable. Form roots before dividing by the total's root.
    root_weights=sqrt(weights):/sqrt(mass)
    weights=weights:/mass
    center=rdistnn_compensated_sum(z:*weights)
    d=root_weights:*(z:-center)
    spread=max(abs(d))
    if(spread==0) return(0)
    z=d:/spread
    norm=spread*sqrt(rdistnn_compensated_sum(z:*z))
    return(norm*scale)
}

real matrix rdistnn_weighted_signals(real matrix signal, real matrix weights)
{
    real matrix out
    real scalar component,target
    if(rows(signal)!=rows(weights)) {
        errprintf("Weighted signals require aligned rows.\n")
        exit(error(198))
    }
    // Materialize each double product before reduction. Unlike quadcross,
    // quadsum(...,1) propagates input and product-overflow missing values.
    out=J(cols(signal),cols(weights),0)
    for(target=1;target<=cols(weights);target++) {
        for(component=1;component<=cols(signal);component++) {
            out[component,target]=quadsum(signal[,component]:*weights[,target],1)
        }
    }
    return(out)
}

void rdistnn_accumulation_add(real matrix high, real matrix low, real matrix value)
{
    real matrix next_high, error, finite, first
    real colvector ix
    real scalar j
    next_high=high+value
    finite=(next_high:<.)
    first=finite:&(abs(high):>=abs(value))
    error=J(rows(high),cols(high),0)
    for(j=1;j<=cols(high);j++) {
        ix=selectindex(first[,j])
        if(length(ix)) error[ix,j]=(high[ix,j]-next_high[ix,j])+value[ix,j]
        ix=selectindex(finite[,j]:&(!first[,j]))
        if(length(ix)) error[ix,j]=(value[ix,j]-next_high[ix,j])+high[ix,j]
    }
    low=low+error
    for(j=1;j<=cols(high);j++) {
        ix=selectindex(!finite[,j])
        if(length(ix)) low[ix,j]=J(length(ix),1,0)
    }
    high=next_high
}

real matrix rdistnn_accumulation_merge(real matrix left_high,
    real matrix right_high, real matrix left_low, real matrix right_low,
    real rowvector offset)
{
    real matrix high,low,value
    high=low=J(rows(left_high),cols(left_high),0)
    rdistnn_accumulation_add(high,low,left_high)
    rdistnn_accumulation_add(high,low,right_high)
    rdistnn_accumulation_add(high,low,left_low)
    rdistnn_accumulation_add(high,low,right_low)
    value=J(rows(high),1,1)*offset
    rdistnn_accumulation_add(high,low,value)
    return(high+low)
}

real matrix rdistnn_accumulation_deletions(real matrix z,
    real colvector table, real rowvector offset)
{
    real scalar n,q,r,last,first,at,component,leading,correction,value,next_high
    real colvector ix
    real rowvector right_high,right_low
    real matrix high,low,suffix_high,suffix_low
    n=rows(z);q=cols(z)
    high=low=J(n,q,0)
    // Scalar scans avoid allocating tiny matrices/masks for every input row.
    // Keep the same Neumaier update and nonfinite policy as accumulation_add.
    for(component=1;component<=q;component++) {
        leading=correction=0
        for(r=2;r<=n;r++) {
            value=table[r-1]*z[r-1,component]
            next_high=leading+value
            if(next_high<.) {
                if(abs(leading)>=abs(value)) correction=correction+((leading-next_high)+value)
                else correction=correction+((value-next_high)+leading)
            }
            else correction=0
            leading=next_high
            high[r,component]=leading;low[r,component]=correction
        }
    }
    right_high=right_low=J(1,q,0)
    last=n
    while(last>0) {
        first=max((1,last-2047));ix=(first::last)
        suffix_high=suffix_low=J(rows(ix),q,0)
        for(component=1;component<=q;component++) {
            leading=right_high[component];correction=right_low[component]
            for(r=last;r>=first;r--) {
                at=r-first+1
                suffix_high[at,component]=leading;suffix_low[at,component]=correction
                if(r>1) {
                    value=table[r-1]*z[r,component]
                    next_high=leading+value
                    if(next_high<.) {
                        if(abs(leading)>=abs(value)) correction=correction+((leading-next_high)+value)
                        else correction=correction+((value-next_high)+leading)
                    }
                    else correction=0
                    leading=next_high
                }
            }
            right_high[component]=leading;right_low[component]=correction
        }
        high[ix,]=rdistnn_accumulation_merge(high[ix,],suffix_high,
            low[ix,],suffix_low,offset)
        last=first-1
    }
    return(high)
}

// Compensated sums are limited to bounded exceptional-scale calculations.
// Products are rounded before entering this shared accumulation convention.
real scalar rdistnn_scale_sum(real colvector x)
{
    real scalar i,total,correction,next_total,value
    total=0; correction=0
    for(i=1;i<=rows(x);i++) {
        value=x[i]; next_total=total+value
        if(abs(total)>=abs(value)) correction=correction+((total-next_total)+value)
        else correction=correction+((value-next_total)+total)
        total=next_total
    }
    return(total+correction)
}

// Mata's ^ uses a transcendental path; construct the binary scale with exact
// multiplications instead. At most eleven iterations cover all double inputs.
real scalar rdistnn_power2(real scalar exponent)
{
    real scalar power,factor,remaining
    factor=exponent<0 ? .5 : 2
    // Mata scalars are passed by reference. Preserve the caller's exponent,
    // including step sizes used by bounded restoration loops.
    remaining=abs(exponent); power=1
    while(remaining>0) {
        if(mod(remaining,2)) power=power*factor
        remaining=floor(remaining/2)
        if(remaining>0) factor=factor*factor
    }
    return(power)
}

real scalar rdistnn_binary_scale(real scalar amplitude)
{
    return(rdistnn_power2(max((-1074,min((1022,floor(ln(amplitude)/ln(2))))))))
}

// Exact-power decomposition shared by bounded range recoveries. The positive
// finite ceiling in Stata is smaller than the permitted negative magnitude.
real rowvector rdistnn_binary_parts(real scalar value)
{
    real scalar amplitude,half,power,normalized
    if(missing(value)) exit(error(198))
    if(value==0) return((0,0))
    amplitude=abs(value); half=0
    if(amplitude>=.) {
        amplitude=abs(value/2); half=1
    }
    power=max((-1074,min((1022,floor(ln(amplitude)/ln(2))+half))))
    normalized=value/rdistnn_power2(power)
    while(abs(normalized)<1 & power>(-1074)) {
        power=power-1
        normalized=value/rdistnn_power2(power)
    }
    while(abs(normalized)>=2 & power<1022) {
        power=power+1
        normalized=value/rdistnn_power2(power)
    }
    return((normalized,power))
}

// Left-continuous weighted empirical quantile used by R diagnostics.
real scalar rdistnn_weighted_quantile(real colvector x,
                                       real colvector w,
                                       real scalar prob)
{
    real scalar index
    real colvector ord

    if (rows(x) != rows(w) | rows(x) < 1) {
        errprintf("rdistnn_weighted_quantile: x and w must have the same positive length.\n")
        exit(error(198))
    }
    if (any(missing(x))) {
        errprintf("rdistnn_weighted_quantile: x must contain only finite, nonmissing values.\n")
        exit(error(198))
    }
    if (missing(prob) | prob < 0 | prob > 1) {
        errprintf("rdistnn_weighted_quantile: probability must be finite and lie in [0,1].\n")
        exit(error(198))
    }
    if (any(missing(w)) | any(w :< 0)) {
        errprintf("rdistnn_weighted_quantile: weights must be finite and nonnegative.\n")
        exit(error(198))
    }
    ord = order((x, (1::rows(x))), (1, 2))
    index=rdistnn_weighted_cdf_index(w[ord],prob)
    if(missing(index)) {
        errprintf("rdistnn_weighted_quantile: total weight must be finite and strictly positive.\n")
        exit(error(198))
    }
    return(x[ord[index]])
}
end
