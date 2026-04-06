*! version 0.0.1  2026-04-05  Jakob Juergens
*! Nonparametric regression at a point via distributional nearest neighbors

program define np_dnn, eclass
    version 16

    if replay() {
        if "`e(cmd)'" != "np_dnn" error 301
        Display `0'
        exit
    }

    syntax varname(numeric) [if] [in] , ///
        x(varlist numeric)               ///  covariates
        x0(numlist)                      ///  evaluation point
        [                                ///
        k(integer 5)                     ///  # nearest neighbors
        s(integer 0)                     ///  kernel order (0 = auto)
        Distance(string)                 ///  l2 (default) | mahalanobis
        Method(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///
        ]

    marksample touse
    markout `touse' `x'

    local outcome `varlist'
    local xvars   `x'
    local kx : word count `xvars'

    quietly count if `touse'
    if r(N) == 0 { di as error "No observations."; exit 2000 }
    local n = r(N)

    local nx0 : word count `x0'
    if `nx0' != `kx' {
        di as error "x0() must have `kx' elements."
        exit 198
    }

    if "`distance'" == "" local distance "l2"
    if "`method'"   == "" local method   "exact"
    if !inlist("`method'", "exact", "asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }
    if `s' == 0  local s = floor(`n'^((`kx' + 1) / (`kx' + 2)))
    if `s' < `k' {
        di as error "s(`s') < k(`k')."
        exit 198
    }

    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    local x0_mata ""
    foreach v of local x0 {
        local x0_mata `x0_mata' `v',
    }
    local x0_mata = substr("`x0_mata'", 1, length("`x0_mata'") - 1)

    // ── plugin fast path for weight computation ───────────────────────
    local __plugin_ok 0
    if ${rdistnn_plugin_ok} {
        tempvar __wvar
        quietly gen double `__wvar' = .
        capture noisily plugin call rdistnn_plugin `xvars' `__wvar' ///
            if `touse', weights_l2 `method' `k' `s' `x0'
        local __plugin_ok = (_rc == 0)
    }

    // ── compute weights ───────────────────────────────────────────────
    if `__plugin_ok'  mata: _np_w = st_data(., "`__wvar'", "`touse'")
    if !`__plugin_ok' mata: _np_w = rdistnn_weights(st_data(., "`xvars'", "`touse'"), (`x0_mata'), `k', `s', "`method'")

    // ── point estimate and LOO jackknife SE ───────────────────────────
    mata: _np_y    = st_data(., "`outcome'", "`touse'")
    mata: _np_mu   = sum(_np_w :* _np_y)
    mata: _np_nn   = rows(_np_w)
    mata: _np_den  = 1 :- _np_w
    mata: _np_safe = (_np_den :> 1e-10)
    mata: _np_loo  = J(_np_nn, 1, _np_mu)
    mata: _np_idx  = selectindex(_np_safe)
    mata: _np_loo[_np_idx] = (_np_mu :- _np_w[_np_idx] :* _np_y[_np_idx]) :/ _np_den[_np_idx]
    mata: _np_sejk = sqrt((_np_nn-1) / _np_nn * sum((_np_loo :- _np_mu):^2))
    mata: _np_neff = 1 / sum(_np_w:^2)

    mata: st_numscalar("r(mu_hat)", _np_mu)
    mata: st_numscalar("r(se_jk)",  _np_sejk)
    mata: st_numscalar("r(neff)",   _np_neff)

    tempname b V
    matrix `b' = (r(mu_hat))
    matrix `V' = (r(se_jk)^2)
    matrix colnames `b' = mu
    matrix rownames `V' = mu
    matrix colnames `V' = mu

    ereturn post `b' `V', esample(`touse')
    ereturn scalar N = `n'
    ereturn scalar neff   = r(neff)
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar kx     = `kx'
    ereturn local  x0     "`x0'"
    ereturn local  xvars  "`xvars'"
    ereturn local  method "`method'"
    ereturn local  cmd    "np_dnn"
    ereturn local  cmdline "np_dnn `0'"

    Display, level(`level')
end

program define Display
    syntax [, Level(cilevel)]
    di ""
    di as text "DNN Nonparametric Regression"
    di as text "Evaluation point x0: " as result "`e(x0)'"
    di as text "k = " as result e(k) as text "  s = " as result e(s) ///
       as text "  neff = " as result %6.2f e(neff)
    di ""
    ereturn display, level(`level')
end
