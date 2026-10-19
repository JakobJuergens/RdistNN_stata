version 17
// Rank laws and observation-level weights; exact tables mirror R and C++.
capture mata: mata drop rdistnn_exact_component()
capture mata: mata drop rdistnn_exact_table()
capture mata: mata drop rdistnn_log_binomial()
capture mata: mata drop rdistnn_weights_exact()
capture mata: mata drop rdistnn_beta_interval()
capture mata: mata drop rdistnn_weights_asymp()
capture mata: mata drop rdistnn_weights_from_ranks()
capture mata: mata drop rdistnn_weights_from_rmat()
capture mata: mata drop rdistnn_weights()
mata:
// ── rdistnn_log_binomial ────────────────────────────────────────────────────────────
// Log binomial coefficient ln(C(n,k)) via lnfactorial().
// Returns -1e300 (effectively -inf) for out-of-support inputs.
real scalar rdistnn_log_binomial(real scalar n, real scalar k)
{
    if (k < 0 | k > n | n < 0) return(-1e300)
    return(lnfactorial(n) - lnfactorial(k) - lnfactorial(n-k))
}

real colvector rdistnn_exact_component(real scalar n, real scalar s, real scalar l)
{
    real scalar last,mode,direction,mantissa,exponent,rank,finish,ratio
    real colvector p
    p=J(n-s+1,1,0); last=n-s+l
    mode=max((l,min((last,floor((n*(l-1))/(s-1))+1))))
    p[mode-l+1]=1
    for(direction=1;direction>=-1;direction=direction-2) {
        mantissa=.5; exponent=1; rank=mode
        if(direction==1) finish=last
        else finish=l
        while(rank!=finish) {
            if(direction==1) ratio=(rank/(rank-l+1))*((n-rank-s+l)/(n-rank))
            else ratio=((rank-l)/(rank-1))*((n-rank+1)/(n-rank-s+l+1))
            mantissa=mantissa*ratio
            while(mantissa<.5) {
                mantissa=mantissa*2
                exponent=exponent-1
            }
            while(mantissa>=1) {
                mantissa=mantissa/2
                exponent=exponent+1
            }
            rank=rank+direction
            if(exponent < -1074) p[rank-l+1]=0
            else p[rank-l+1]=mantissa*rdistnn_power2(exponent)
        }
    }
    return(p/rdistnn_compensated_sum(p))
}

real colvector rdistnn_exact_table(real scalar n, real scalar k, real scalar s)
{
    real scalar l
    real colvector high,low,p,out,at,h,b
    if(k==s) return(J(n,1,1/n))
    high=low=J(n,1,0)
    for(l=1;l<=k;l++) {
        p=rdistnn_exact_component(n,s,l)
        at=(l::(n-s+l))
        h=high[at]; b=low[at]
        rdistnn_accumulation_add(h,b,p)
        high[at]=h; low[at]=b
    }
    out=(high+low)/k
    return(out/rdistnn_compensated_sum(out))
}

// ── rdistnn_weights_exact ─────────────────────────────────────────────
// Exact DNN rank weights (Steele closed form / full WOR hypergeometric).
//
// For obs with global rank r (1 = nearest to x0):
//
//   w_r = (1/k) * sum_{l=1}^{k} C(r-1, l-1) * C(n-r, s-l) / C(n, s)
//
// Positive-weight support: r = 1, ..., n-s+k  (NOT 1..s).  Thus larger s
// shrinks support and tightens localization; for k=1 the edge is n-s+1.
// Here k is the number selected within each s-subset, not the global support
// size.  For k=1 the weight equals C(n-r, s-1) / C(n, s), the Steele formula.
//
// Arguments:
//   ranks : n-vector of integer ranks (1 = nearest); n = rows(ranks)
//   k     : number selected within each s-subset (integer >= 1)
//   s     : kernel order / subset size (integer, k <= s <= n)
//
// Returns: n-vector of weights summing to 1 (up to floating-point error).
//
real colvector rdistnn_weights_exact(real colvector ranks,
                                      real scalar k,
                                      real scalar s)
{
    real colvector table
    real scalar n
    n=rows(ranks)
    if(k<1 | k>s | s>n | k!=floor(k) | s!=floor(s)) {
        errprintf("rdistnn_weights_exact: require 1 <= k <= s <= n.\n")
        exit(error(198))
    }
    table=rdistnn_exact_table(n,k,s)
    return(table[ranks])
}

// Stable Beta interval probability; use the small upper tail directly above
// the shared lower-endpoint crossover. Preserve every representable mass.
real scalar rdistnn_beta_interval(real scalar a, real scalar b,
    real scalar lower, real scalar upper)
{
    if(lower<=a/(a+b)) return(ibeta(a,b,upper)-ibeta(a,b,lower))
    return(ibetatail(a,b,lower)-ibetatail(a,b,upper))
}

// ── rdistnn_weights_asymp ─────────────────────────────────────────────
// Asymptotic DNN weights via Beta order-statistic increments.
//
// For rank j (1-based), the asymptotic weight is:
//   w_rank[j] = (1/k) * sum_{l=1}^{k} [ F_Beta(l, s-l+1)(j/n)
//                                       - F_Beta(l, s-l+1)((j-1)/n) ]
//
// This is the large-n approximation to the exact WOR formula, treating
// the sorted s-subset positions as i.i.d. Uniform(0,1) order statistics.
// Full support over 1..n; renormalized to sum exactly 1.
//
// Arguments:
//   ranks : n-vector of integer ranks (1 = nearest)
//   k     : number selected within each s-subset (integer >= 1)
//   s     : kernel order (integer >= k)
//
// Returns: n-vector of weights summing to 1.
//
real colvector rdistnn_weights_asymp(real colvector ranks,
                                      real scalar k,
                                      real scalar s)
{
    real scalar     n, i, l, r
    real colvector  w, w_rank

    n      = rows(ranks)
    w_rank = J(n, 1, 0)

    for (l = 1; l <= k; l++) {
        // Beta(l, s-l+1) CDF increments on ((j-1)/n, j/n] for j=1..n
        for (i = 1; i <= n; i++) {
            w_rank[i] = w_rank[i] + rdistnn_beta_interval(l, s-l+1, (i-1)/n, i/n)
        }
    }
    w_rank = w_rank / k

    // Normalize for safety
    w_rank = w_rank / rdistnn_sum(w_rank)

    // Map from rank order back to observation order
    w = J(n, 1, 0)
    for (i = 1; i <= n; i++) {
        r    = ranks[i]
        w[i] = w_rank[r]
    }
    return(w)
}

// ── rdistnn_weights_from_ranks ────────────────────────────────────────
// Compute DNN weights directly from a precomputed rank vector.
// Factored out of rdistnn_weights() so the jackknife can reuse it
// without recomputing distances.
//
// Arguments:
//   ranks  : n-vector of integer ranks (1 = nearest)
//   k      : number selected within each s-subset
//   s      : kernel order
//   method : "exact" (default) or "asymp"
//
// Returns: n-vector of DNN weights.
//
real colvector rdistnn_weights_from_ranks(real colvector ranks,
                                           real scalar k,
                                           real scalar s,
                                           | string scalar method)
{
    if (args() < 4 | method == "" | method == "exact") {
        return(rdistnn_weights_exact(ranks, k, s))
    }
    else if (method == "asymp") {
        return(rdistnn_weights_asymp(ranks, k, s))
    }
    else {
        errprintf("rdistnn_weights_from_ranks: unknown method '%s'. Use 'exact' or 'asymp'.\n", method)
        exit(error(198))
    }
}

// Apply the rank-to-weight map independently to every target/rank column.
real matrix rdistnn_weights_from_rmat(real matrix ranks,
                                       real scalar k,
                                       real scalar s,
                                       | string scalar method)
{
    real scalar j, m
    real colvector table
    real matrix W

    if (args() < 4 | method == "") method = "exact"
    m = cols(ranks)
    if (m < 1) {
        errprintf("rdistnn_weights_from_rmat: ranks must have at least one column.\n")
        exit(error(198))
    }

    // The rank-weight table is independent of the target permutation.
    W = J(rows(ranks), m, .)
    table = rdistnn_weights_from_ranks((1::rows(ranks)), k, s, method)
    for (j = 1; j <= m; j++) {
        W[,j] = table[ranks[,j]]
    }
    return(W)
}

// ── rdistnn_weights ───────────────────────────────────────────────────
// Top-level: compute DNN weights for X evaluated at x0.
//
// Arguments:
//   X        : n x kx matrix
//   x0       : 1 x kx evaluation point
//   k        : number selected within each s-subset
//   s        : kernel order
//   method   : "exact" (default) or "asymp"
//   distance : "l2" (default) or "mahalanobis"
//   Sinv     : precomputed inverse covariance (optional; computed if not supplied
//              and distance == "mahalanobis")
//
// Returns: n-vector of DNN weights.
//
real colvector rdistnn_weights(real matrix X, real rowvector x0,
                                real scalar k, real scalar s,
                                | string scalar method,
                                  string scalar distance,
                                  real matrix Sinv)
{
    real colvector dist, ranks

    if (args() < 5 | method == "")   method   = "exact"
    if (args() < 6 | distance == "") distance = "l2"

    if (args() < 7 | rows(Sinv) == 0) Sinv=rdistnn_fit_geometry(X,distance,"sd",J(1,0,.),1e-12)
    dist=rdistnn_distance_matrix(X,x0,distance,Sinv)
    ranks = rdistnn_order(dist)
    return(rdistnn_weights_from_ranks(ranks, k, s, method))
}
end
