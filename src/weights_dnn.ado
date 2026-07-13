*! version 0.0.7  2026-07-12  Jakob Juergens
*! DNN rank weights at a fixed covariate point x0

program define weights_dnn, rclass
    version 17

    syntax varlist(min=1) [if] [in] ,   ///
        x0(numlist)                      ///  evaluation point (length = # vars)
        [                                ///
        k(integer 1)                     ///  # selected within each s-subset
        s(integer 0)                     ///  subset/kernel order; larger = tighter
        Generate(name)                   ///  new variable to store weights
        Distance(string)                 ///  l2 (default) | mahalanobis
        Method(string)                   ///  exact (default) | asymp
        ]

    // ── mark sample ──────────────────────────────────────────────────
    marksample touse
    markout `touse' `varlist'

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

    if "`distance'" == "" local distance "l2"
    if !inlist("`distance'", "l2", "mahalanobis") {
        di as error "distance() must be l2 or mahalanobis."
        exit 198
    }

    // Mahalanobis requires at least 2 covariates to form a covariance matrix
    if "`distance'" == "mahalanobis" & `kx' < 2 {
        di as error "distance(mahalanobis) requires at least 2 covariates."
        exit 198
    }

    if "`method'" == "" local method "exact"
    if !inlist("`method'", "exact", "asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }

    // Paper finite-sample default: floor(0.4*n^((kx+1)/(kx+2)))
    quietly count if `touse'
    local n = r(N)
    if `s' == 0 {
        local s = floor(.4 * `n'^((`kx' + 1) / (`kx' + 2)))
    }
    if `s' < `k' {
        di as error "s(`s') < k(`k'); increase s or decrease k."
        exit 198
    }
    if `s' > `n' {
        di as error "s(`s') must be no greater than n(`n')."
        exit 198
    }

    // ── load Mata library (also attempts plugin declaration) ──────────
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    // ── build x0 string for Mata ──────────────────────────────────────
    local x0_mata ""
    foreach v of local x0 {
        local x0_mata `x0_mata' `v',
    }
    local x0_mata = substr("`x0_mata'", 1, length("`x0_mata'") - 1)

    // ── build plugin-safe x0 arglist ─────────────────────────────────
    // Negative numeric args can be split at the plugin-call arglist
    // boundary unless each component is passed as one quoted string.
    local dq = char(34)
    local x0_plugin ""
    foreach v of local x0 {
        local x0_plugin `"`x0_plugin' `dq'`v'`dq'"'
    }

    // ── compute weights ───────────────────────────────────────────────
    // Eligible L2 calls try the installed plugin on every platform. Any load
    // or execution failure falls back to Mata unless the plugin is required.
    local __plugin_ok 0
    local __plugin_mode "${rdistnn_use_plugin}"
    local __require_plugin = ("`__plugin_mode'" == "1")
    local __try_plugin = ("`distance'" == "l2" & "`__plugin_mode'" != "0")

    if `__try_plugin' & !${rdistnn_plugin_ok} & `__require_plugin' {
        di as error "rdistnn.plugin was requested but is not available."
        exit 601
    }

    if `__try_plugin' & ${rdistnn_plugin_ok} {
        tempvar __wvar
        quietly gen double `__wvar' = .
        if `__require_plugin' {
            capture noisily plugin call rdistnn_plugin `xvars' `__wvar' ///
                if `touse', weights_l2 `method' `k' `s' `x0_plugin'
        }
        else {
            capture plugin call rdistnn_plugin `xvars' `__wvar' ///
                if `touse', weights_l2 `method' `k' `s' `x0_plugin'
        }
        local __plugin_rc = _rc
        local __plugin_ok = (`__plugin_rc' == 0)
        if !`__plugin_ok' & `__require_plugin' {
            exit `__plugin_rc'
        }
        if `__plugin_ok' {
            mata: st_local("wt_neff", strofreal(1 / sum(st_data(., "`__wvar'", "`touse'") :^2), "%21.0g"))
            if "`generate'" != "" {
                quietly gen double `generate' = `__wvar' if `touse'
                di as text "DNN weights stored in `generate'."
            }
        }
    }

    if `__plugin_ok' {
        return scalar neff = `wt_neff'
        return scalar n    = `n'
        return scalar k    = `k'
        return scalar s    = `s'
        return scalar kx   = `kx'

        di as text "Effective sample size (neff): " as result %8.2f `wt_neff'
        exit
    }

    // Pure-Mata fallback — handles both L2 and Mahalanobis.
    // rdistnn_ranks() does the L2/Mahal dispatch internally to avoid
    // if-else braces in the inline mata { } block (Stata parser limit).
    mata {
        _X     = st_data(., "`xvars'", "`touse'")
        _x0    = (`x0_mata')
        _ranks = rdistnn_ranks(_X, _x0, "`distance'")
        _w     = rdistnn_weights_from_ranks(_ranks, `k', `s', "`method'")
        st_local("wt_neff", strofreal(1 / sum(_w :^2), "%21.0g"))
    }
    if "`generate'" != "" {
        quietly gen double `generate' = .
        mata: st_store(., "`generate'", "`touse'", _w)
        di as text "DNN weights stored in `generate'."
    }

    return scalar neff = `wt_neff'
    return scalar n    = `n'
    return scalar k    = `k'
    return scalar s    = `s'
    return scalar kx   = `kx'

    di as text "Effective sample size (neff): " as result %8.2f `wt_neff'

end
