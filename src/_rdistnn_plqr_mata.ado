version 17
// Adaptive PLQR predicate: retain certified ordinary comparisons, otherwise
// sum the binary64 products exactly with separate base-2^26 magnitudes.
foreach fn in indicator exact_indicator exact_product exact_piece {
    capture mata: mata drop rdistnn_plqr_`fn'()
}
mata:
void rdistnn_plqr_exact_piece(real colvector limbs, real scalar piece,
    real scalar bit)
{
    real scalar index,value,next_value,carry
    if(piece==0) return
    if(bit<0) exit(error(498))
    index=floor(bit/26)+1
    value=piece*rdistnn_power2(mod(bit,26))
    while(value!=0) {
        if(index>rows(limbs)) {
            errprintf("Internal PLQR exact comparison capacity exceeded.\n")
            exit(error(498))
        }
        next_value=limbs[index]+value
        carry=floor(next_value/67108864)
        limbs[index]=next_value-carry*67108864
        value=carry; index=index+1
    }
}

void rdistnn_plqr_exact_product(real colvector positive,
    real colvector negative, real scalar a, real scalar b)
{
    real rowvector left,right,xs,ys
    real scalar x,y,xnext,ynext,i,j,product,high,low,start,bit
    if(a==0 | b==0) return
    left=rdistnn_binary_parts(a); right=rdistnn_binary_parts(b)
    x=abs(left[1])*4503599627370496
    y=abs(right[1])*4503599627370496
    xs=ys=J(1,3,0)
    for(i=1;i<=3;i++) {
        xnext=floor(x/67108864); ynext=floor(y/67108864)
        xs[i]=x-xnext*67108864; ys[i]=y-ynext*67108864
        x=xnext; y=ynext
    }
    // Origin -2252 covers two minimum operand powers, with the two
    // 52-bit significand shifts. Split products before the bit shift so
    // every multiplication/addition is an exact integer below 2^52.
    start=left[2]+right[2]-104-(-2252)
    for(i=1;i<=3;i++) for(j=1;j<=3;j++) {
        product=xs[i]*ys[j]
        high=floor(product/67108864); low=product-high*67108864
        bit=start+26*(i+j-2)
        if((a<0)==(b<0)) {
            rdistnn_plqr_exact_piece(positive,low,bit)
            rdistnn_plqr_exact_piece(positive,high,bit+26)
        }
        else {
            rdistnn_plqr_exact_piece(negative,low,bit)
            rdistnn_plqr_exact_piece(negative,high,bit+26)
        }
    }
}

real scalar rdistnn_plqr_exact_indicator(real rowvector D,
    real colvector theta, real scalar g, real scalar y)
{
    real colvector positive,negative
    real scalar j
    // 170 limbs cover [-2252,2168). At most INT_MAX+2 products, each
    // smaller than 2^2048, have magnitude below 2^2080 in either sign.
    positive=negative=J(170,1,0)
    for(j=1;j<=cols(D);j++) {
        rdistnn_plqr_exact_product(positive,negative,D[j],theta[j])
    }
    rdistnn_plqr_exact_product(positive,negative,g,1)
    // Keep y and -1 separate: negating a finite negative extreme can
    // exceed Stata's positive ceiling before the exact product is formed.
    rdistnn_plqr_exact_product(positive,negative,y,-1)
    for(j=170;j>=1;j--) {
        if(positive[j]>negative[j]) return(1)
        if(positive[j]<negative[j]) return(0)
    }
    return(1)
}

real colvector rdistnn_plqr_indicator(real matrix D, real colvector theta,
    real colvector g, real colvector y)
{
    real scalar n,q,i,j,risk,product,margin,bound,gate,tiny
    real colvector threshold,below
    n=rows(D); q=cols(D)
    if(n<1 | q<1 | q>2147483647 | rows(theta)!=q | rows(g)!=n |
        rows(y)!=n | hasmissing(D) | hasmissing(theta) |
        hasmissing(g) | hasmissing(y)) {
        errprintf("PLQR comparisons require finite aligned regressors, coefficients, nuisance values and outcomes.\n")
        exit(error(198))
    }
    gate=rdistnn_power2(1022); tiny=rdistnn_power2(-1074)
    threshold=D*theta:+g
    below=(y:<=threshold)
    for(i=1;i<=n;i++) {
        // Match the R/native risk bound and materialized product order.
        risk=abs(g[i])+abs(y[i])
        for(j=1;j<=q;j++) {
            product=D[i,j]*theta[j]
            risk=risk+abs(product)
        }
        if(risk<gate & threshold[i]<.) {
            margin=threshold[i]-y[i]
            bound=8*(q+3)*epsilon(1)*risk+32*(q+2)*tiny
            if(margin<. & abs(margin)>bound) continue
        }
        below[i]=rdistnn_plqr_exact_indicator(D[i,],theta,g[i],y[i])
    }
    return(below)
}
end
