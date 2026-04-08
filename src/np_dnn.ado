*! version 0.0.2  2026-04-06  Jakob Juergens
*! Nonparametric regression at a point via distributional nearest neighbors

program define np_dnn, eclass
    version 16

    if replay() {
        if "`e(cmd)'" != "np_dnn" error 301
        _np_dnn_Display `0'
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
    markout `touse' `x' `varlist'

    local outcome `varlist'
    local xvars   `x'
    local kx : word count `xvars'

    quietly count if `touse'
    if r(N) == 0 {
        di as error "No observations."
        exit 2000
    }
    local n = r(N)

    local nx0 : word count `x0'
    if `nx0' != `kx' {
        di as error "x0() must have `kx' elements."
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

    if `s' == 0  local s = floor(`n'^((`kx' + 1) / (`kx' + 2)))
    if `s' < `k' {
        di as error "s(`s') < k(`k')."
        exit 198
    }
    if `s' >= `n' {
        di as error "s(`s') must be less than n(`n') for jackknife SE."
        exit 198
    }

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
        _np_X     = st_data(., "`xvars'", "`touse'")
        _np_x0    = (`x0_mata')
        _np_ranks = rdistnn_ranks(_np_X, _np_x0, "`distance'")
        _np_w     = rdistnn_weights_from_ranks(_np_ranks, `k', `s', "`method'")
        _np_y     = st_data(., "`outcome'", "`touse'")
        _np_mu    = sum(_np_w :* _np_y)
        _np_sejk  = rdistnn_jk_se_from_ranks(_np_ranks, _np_y, `k', `s', "`method'")
        _np_neff  = 1 / sum(_np_w :^2)
        st_matrix("`b'", _np_mu)
        st_matrix("`V'", _np_sejk^2)
        st_numscalar("r(neff)", _np_neff)
    }
    matrix colnames `b' = mu
    matrix rownames `V' = mu
    matrix colnames `V' = mu

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
    ereturn local  cmd      "np_dnn"
    ereturn local  cmdline  "np_dnn `0'"

    _np_dnn_Display, level(`level')
end

program define _np_dnn_Display
    syntax [, Level(cilevel)]
    di ""
    di as text "DNN Nonparametric Regression"
    di as text "Evaluation point x0: " as result "`e(x0)'"
    di as text "Distance: " as result "`e(distance)'" ///
       as text "  k = " as result e(k) ///
       as text "  s = " as result e(s) ///
       as text "  neff = " as result %6.2f e(neff)
    di ""
    ereturn display, level(`level')
end
