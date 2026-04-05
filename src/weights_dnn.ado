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
        Method(string)                   ///  exact (default) | asymp
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

    if "`method'" == "" local method "exact"
    if !inlist("`method'", "exact", "asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }

    // ── load Mata library (also attempts plugin declaration) ──────────
    quietly run "_rdistnn_mata.ado"

    // ── build x0 string for Mata and plugin ───────────────────────────
    local x0_mata ""
    foreach v of local x0 {
        local x0_mata `x0_mata' `v',
    }
    local x0_mata = substr("`x0_mata'", 1, length("`x0_mata'") - 1)

    // ── compute weights (plugin fast path, then Mata fallback) ────────
    // Weights are written to a temp var (or directly to generate if set).
    // The plugin reads all observations satisfying `touse'; Mata does the
    // same via st_data.

    if ${rdistnn_plugin_ok} {
        // create the output weight variable up front
        if "`generate'" != "" {
            quietly gen double `generate' = .
            capture noisily plugin call rdistnn_plugin `xvars' `generate' ///
                if `touse', weights_l2 `method' `k' `s' `x0'
        }
        else {
            tempvar __wvar
            quietly gen double `__wvar' = .
            capture noisily plugin call rdistnn_plugin `xvars' `__wvar' ///
                if `touse', weights_l2 `method' `k' `s' `x0'
        }
        local __plugin_ok = (_rc == 0)
    }
    else {
        local __plugin_ok 0
    }

    if `__plugin_ok' {
        // neff from the stored weight variable
        if "`generate'" != "" {
            mata: st_numscalar("r(neff)", 1 / sum(st_data(., "`generate'", "`touse'") :^2))
            di as text "DNN weights stored in `generate'."
        }
        else {
            mata: st_numscalar("r(neff)", 1 / sum(st_data(., "`__wvar'", "`touse'") :^2))
        }
    }
    else {
        // Pure-Mata fallback
        mata {
            _X = st_data(., "`xvars'", "`touse'")
            _x0 = (`x0_mata')
            _w  = rdistnn_weights(_X, _x0, `k', `s', "`method'")
        }
        if "`generate'" != "" {
            quietly gen double `generate' = .
            mata: st_store(., "`generate'", "`touse'", _w)
            di as text "DNN weights stored in `generate'."
        }
        mata: st_numscalar("r(neff)", 1 / sum(_w :^2))
    }
    return scalar neff = r(neff)
    return scalar n    = `n'
    return scalar k    = `k'
    return scalar s    = `s'
    return scalar kx   = `kx'

    di as text "Effective sample size (neff): " as result %8.2f r(neff)

end
