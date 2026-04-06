*! version 0.0.1  2026-04-05  Jakob Juergens
*! Conditional average treatment effect via distributional nearest neighbors

program define cate_dnn, eclass
    version 16

    // ── replay ────────────────────────────────────────────────────────
    if replay() {
        if "`e(cmd)'" != "cate_dnn" error 301
        Display `0'
        exit
    }

    syntax varname(numeric) [if] [in] , ///
        TReat(varname numeric)           ///  binary treatment D
        x(varlist numeric)               ///  covariates
        x0(numlist)                      ///  evaluation point
        [                                ///
        Propensity(varname numeric)      ///  estimated/oracle e(X)
        M0(varname numeric)              ///  estimated/oracle E[Y|X,D=0]
        M1(varname numeric)              ///  estimated/oracle E[Y|X,D=1]
        k(integer 5)                     ///  # nearest neighbors
        s(integer 0)                     ///  kernel order (0 = auto)
        Distance(string)                 ///  l2 (default) | mahalanobis
        Method(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///  CI level (default 95)
        noCONStant                       ///  (reserved; absorbed silently)
        ]

    // ── mark sample ──────────────────────────────────────────────────
    marksample touse
    markout `touse' `treat' `x' `propensity' `m0' `m1'

    local outcome `varlist'
    local xvars   `x'
    local kx : word count `xvars'

    // ── validate ──────────────────────────────────────────────────────
    quietly count if `touse'
    if r(N) == 0 {
        di as error "No observations."
        exit 2000
    }
    local n = r(N)

    local nx0 : word count `x0'
    if `nx0' != `kx' {
        di as error "x0() must have `kx' elements (one per covariate in x())."
        exit 198
    }

    if `k' < 1 {
        di as error "k() must be >= 1."
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

    // ── default s ────────────────────────────────────────────────────
    if `s' == 0 {
        local s = floor(`n'^((`kx' + 1) / (`kx' + 2)))
    }
    if `s' < `k' {
        di as error "s(`s') < k(`k'). Increase n, increase s(), or decrease k()."
        exit 198
    }

    // ── pseudo-outcome ────────────────────────────────────────────────
    // Robinson (1988) pseudo-outcome:  rho_i = (D_i - e_i)(Y_i - m_i) / (e_i(1-e_i))
    // Requires propensity + one of {m0, m1} or their combination m = m0 + e*(m1-m0).
    if "`propensity'" == "" {
        di as error "propensity() is required.  Provide oracle or estimated propensity scores."
        exit 198
    }
    if "`m0'" == "" & "`m1'" == "" {
        di as error "m0() or m1() (or both) is required."
        exit 198
    }

    // Construct m = E[Y|X] = m0 + e*(m1-m0)
    tempvar e_var m_var rho_var w_var

    quietly gen double `e_var' = `propensity' if `touse'

    if "`m0'" != "" & "`m1'" != "" {
        quietly gen double `m_var' = `m0' + `e_var' * (`m1' - `m0') if `touse'
    }
    else if "`m0'" != "" {
        // m = m0 + e*(Y_bar - m0) approximation -- user should supply m1
        di as text "(note: m1() not supplied; using m_i = m0_i as outcome model)"
        quietly gen double `m_var' = `m0' if `touse'
    }
    else {
        di as text "(note: m0() not supplied; using m_i = m1_i as outcome model)"
        quietly gen double `m_var' = `m1' if `touse'
    }

    quietly gen double `rho_var' = ///
        (`treat' - `e_var') * (`outcome' - `m_var') / ///
        max(`e_var' * (1 - `e_var'), 1e-8) if `touse'

    // ── load Mata and compute weights ─────────────────────────────────
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
    if `__plugin_ok'  mata: _ct_w = st_data(., "`__wvar'", "`touse'")
    if !`__plugin_ok' mata: _ct_w = rdistnn_weights(st_data(., "`xvars'", "`touse'"), (`x0_mata'), `k', `s', "`method'")

    // ── point estimate and LOO jackknife SE ───────────────────────────
    mata: _ct_rho  = st_data(., "`rho_var'", "`touse'")
    mata: _ct_tau  = sum(_ct_w :* _ct_rho)
    mata: _ct_nn   = rows(_ct_w)
    mata: _ct_den  = 1 :- _ct_w
    mata: _ct_safe = (_ct_den :> 1e-10)
    mata: _ct_loo  = J(_ct_nn, 1, _ct_tau)
    mata: _ct_idx  = selectindex(_ct_safe)
    mata: _ct_loo[_ct_idx] = (_ct_tau :- _ct_w[_ct_idx] :* _ct_rho[_ct_idx]) :/ _ct_den[_ct_idx]
    mata: _ct_sejk = sqrt((_ct_nn-1) / _ct_nn * sum((_ct_loo :- _ct_tau):^2))
    mata: _ct_neff = 1 / sum(_ct_w:^2)

    mata: st_numscalar("r(tau_hat)", _ct_tau)
    mata: st_numscalar("r(se_jk)",   _ct_sejk)
    mata: st_numscalar("r(neff)",    _ct_neff)

    // ── store eclass results ─────────────────────────────────────────
    tempname b V
    matrix `b' = (r(tau_hat))
    matrix `V' = (r(se_jk)^2)
    matrix colnames `b' = tau
    matrix rownames `V' = tau
    matrix colnames `V' = tau

    ereturn post `b' `V', esample(`touse')
    ereturn scalar N = `n'
    ereturn scalar neff   = r(neff)
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar kx     = `kx'
    ereturn local  x0       "`x0'"
    ereturn local  xvars    "`xvars'"
    ereturn local  distance "`distance'"
    ereturn local  method   "`method'"
    ereturn local  cmd    "cate_dnn"
    ereturn local  cmdline "cate_dnn `0'"

    Display, level(`level')
end

// ── display subroutine ────────────────────────────────────────────────
program define Display
    syntax [, Level(cilevel)]

    di ""
    di as text "DNN Conditional Average Treatment Effect Estimator"
    di as text "Evaluation point x0: " as result "`e(x0)'"
    di as text "Covariates: "          as result "`e(xvars)'"
    di as text "k = " as result e(k) as text "  s = " as result e(s) ///
       as text "  neff = " as result %6.2f e(neff)
    di ""
    ereturn display, level(`level')
end
