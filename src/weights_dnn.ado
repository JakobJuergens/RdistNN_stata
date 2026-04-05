*! version 0.0.1  2026-04-05  Jakob Juergens
*! DNN rank weights at a fixed covariate point x0

program define weights_dnn, rclass
    version 16

    syntax varlist(min=2) [if] [in] ,   ///
        x0(numlist)                      ///  evaluation point (length = # vars)
        k(integer)                       ///  # nearest neighbors
        [                                ///
        s(integer 0)                     ///  kernel order (0 = auto)
        Generate(name)                   ///  new variable to store weights
        Distance(string)                 ///  l2 (default) | mahalanobis
        ]

    // ── mark sample ──────────────────────────────────────────────────
    marksample touse

    // ── parse and validate ────────────────────────────────────────────
    local xvars `varlist'
    local kx    : word count `xvars'

    if `k' < 1 {
        di as error "k() must be a positive integer."
        exit 198
    }

    local nx0 : word count `x0'
    if `nx0' != `kx' {
        di as error "x0() must have `kx' elements (one per covariate), " ///
                    "got `nx0'."
        exit 198
    }

    // Default s: floor(n^((kx+1)/(kx+2)))
    quietly count if `touse'
    local n = r(N)
    if `s' == 0 {
        local s = floor(`n'^((`kx' + 1) / (`kx' + 2)))
    }
    if `s' < `k' {
        di as error "s(`s') < k(`k'); increase s or decrease k."
        exit 198
    }

    if "`distance'" == "" local distance "l2"
    if !inlist("`distance'", "l2", "mahalanobis") {
        di as error "distance() must be l2 or mahalanobis."
        exit 198
    }

    // ── load Mata library ─────────────────────────────────────────────
    quietly run "_rdistnn_mata.ado"

    // ── build covariate matrix in Mata ────────────────────────────────
    // Parse x0 numlist into a Mata row vector
    local x0_mata ""
    foreach v of local x0 {
        local x0_mata `x0_mata' `v',
    }
    local x0_mata = substr("`x0_mata'", 1, length("`x0_mata'") - 1)

    mata {
        _X   = st_data(., "`xvars'", "`touse'")
        _x0  = (`x0_mata')
        _k   = `k'
        _s   = `s'
        _w   = rdistnn_weights(_X, _x0, _k, _s)
    }

    // ── store weights ─────────────────────────────────────────────────
    if "`generate'" != "" {
        quietly gen double `generate' = .
        mata: st_store(., "`generate'", "`touse'", _w)
        di as text "DNN weights stored in `generate'."
    }

    // ── return scalars ────────────────────────────────────────────────
    mata: st_numscalar("r(neff)", 1 / sum(_w :^2))
    return scalar neff = r(neff)
    return scalar n    = `n'
    return scalar k    = `k'
    return scalar s    = `s'
    return scalar kx   = `kx'

    di as text "Effective sample size (neff): " as result %8.2f r(neff)

end
