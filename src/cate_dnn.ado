*! version 0.0.2  2026-04-06  Jakob Juergens
*! Conditional average treatment effect via distributional nearest neighbors

program define cate_dnn, eclass
    version 16

    // ── replay ────────────────────────────────────────────────────────
    if replay() {
        if "`e(cmd)'" != "cate_dnn" error 301
        _cate_dnn_Display `0'
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
    markout `touse' `varlist' `treat' `x' `propensity' `m0' `m1'

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
    if `s' >= `n' {
        di as error "s(`s') must be less than n(`n') for jackknife SE."
        exit 198
    }

    // ── pseudo-outcome ────────────────────────────────────────────────
    // Full AIPW (doubly robust) when both m0 and m1 are provided:
    //   phi_i = D*(Y-m1)/den - (1-D)*(Y-m0)/den + (m1-m0)
    // where den = max(e*(1-e), 1e-8).
    //
    // Fallback to Robinson form when only one nuisance is available:
    //   phi_i = (D-e)*(Y-m)/den
    if "`propensity'" == "" {
        di as error "propensity() is required.  Provide oracle or estimated propensity scores."
        exit 198
    }
    if "`m0'" == "" & "`m1'" == "" {
        di as error "m0() or m1() (or both) is required."
        exit 198
    }

    tempvar den_var rho_var

    quietly gen double `den_var' = max(`propensity' * (1 - `propensity'), 1e-8) if `touse'

    if "`m0'" != "" & "`m1'" != "" {
        // Full AIPW pseudo-outcome
        quietly gen double `rho_var' = ///
            `treat' * (`outcome' - `m1') / `den_var'       ///
            - (1 - `treat') * (`outcome' - `m0') / `den_var' ///
            + (`m1' - `m0')                                ///
            if `touse'
    }
    else {
        // Single nuisance: Robinson form  phi_i = (D-e)*(Y-m)/den
        if "`m0'" != "" {
            di as text "(note: m1() not supplied; using Robinson form with m0)"
            local m_use `m0'
        }
        else {
            di as text "(note: m0() not supplied; using Robinson form with m1)"
            local m_use `m1'
        }
        quietly gen double `rho_var' = ///
            (`treat' - `propensity') * (`outcome' - `m_use') / `den_var' ///
            if `touse'
    }

    // ── load Mata ────────────────────────────────────────────────────
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    local x0_mata ""
    foreach v of local x0 {
        local x0_mata `x0_mata' `v',
    }
    local x0_mata = substr("`x0_mata'", 1, length("`x0_mata'") - 1)

    // ── compute ranks, weights, point estimate, and true delete-1 JK SE ─
    // Always use the Mata path so ranks are available for the JK.
    // rdistnn_ranks() handles L2/Mahalanobis dispatch internally, avoiding
    // if-else braces inside the inline mata { } block (Stata parser limit).
    tempname b V
    mata {
        _ct_X     = st_data(., "`xvars'", "`touse'")
        _ct_x0    = (`x0_mata')
        _ct_ranks = rdistnn_ranks(_ct_X, _ct_x0, "`distance'")
        _ct_w     = rdistnn_weights_from_ranks(_ct_ranks, `k', `s', "`method'")
        _ct_rho   = st_data(., "`rho_var'", "`touse'")
        _ct_tau   = sum(_ct_w :* _ct_rho)
        _ct_sejk  = rdistnn_jk_se_from_ranks(_ct_ranks, _ct_rho, `k', `s', "`method'")
        _ct_neff  = 1 / sum(_ct_w :^2)
        st_matrix("`b'", _ct_tau)
        st_matrix("`V'", _ct_sejk^2)
        st_numscalar("r(neff)", _ct_neff)
    }

    // ── store eclass results ─────────────────────────────────────────
    matrix colnames `b' = tau
    matrix rownames `V' = tau
    matrix colnames `V' = tau

    ereturn post `b' `V', esample(`touse')
    ereturn scalar N      = `n'
    ereturn scalar neff   = scalar(r(neff))
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar kx     = `kx'
    ereturn local  x0       "`x0'"
    ereturn local  xvars    "`xvars'"
    ereturn local  distance "`distance'"
    ereturn local  method   "`method'"
    ereturn local  cmd      "cate_dnn"
    ereturn local  cmdline  "cate_dnn `0'"

    _cate_dnn_Display, level(`level')
end

// ── display subroutine ────────────────────────────────────────────────
program define _cate_dnn_Display
    syntax [, Level(cilevel)]

    di ""
    di as text "DNN Conditional Average Treatment Effect Estimator"
    di as text "Evaluation point x0: " as result "`e(x0)'"
    di as text "Covariates: "          as result "`e(xvars)'"
    di as text "Distance: " as result "`e(distance)'" ///
       as text "  k = " as result e(k) ///
       as text "  s = " as result e(s) ///
       as text "  neff = " as result %6.2f e(neff)
    di ""
    ereturn display, level(`level')
end
