*! version 0.0.1  2026-05-30  Jakob Juergens
*! PLQR orthogonal moments at a point via distributional nearest neighbors

program define plqr_dnn, eclass
    version 17

    if replay() {
        if "`e(cmd)'" != "plqr_dnn" error 301
        _plqr_dnn_Display `0'
        exit
    }

    syntax varname(numeric) [if] [in] , ///
        x(varlist numeric)               ///  covariates
        x0(numlist)                      ///  evaluation point
        D(varlist numeric)               ///  partially linear regressors
        Q(varname numeric)               ///  supplied q_tau(X)
        MD(varlist numeric)              ///  supplied m_D(X), same order as d()
        beta(numlist)                    ///  fixed beta vector
        [                                ///
        tau(real 0.5)                    ///  quantile index
        k(integer 5)                     ///  # nearest neighbors
        s(integer 0)                     ///  kernel order (0 = auto)
        Distance(string)                 ///  l2 (default) | mahalanobis
        Method(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///
        ]

    marksample touse
    markout `touse' `varlist' `x' `d' `q' `md'

    local outcome `varlist'
    local xvars   `x'
    local dvars   `d'
    local mdvars  `md'
    local qvar    `q'
    local kx : word count `xvars'
    local kd : word count `dvars'
    local kmd : word count `mdvars'
    local nbeta : word count `beta'

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

    if `kd' != `kmd' {
        di as error "d() and md() must contain the same number of variables."
        exit 198
    }

    if `nbeta' != `kd' {
        di as error "beta() must have `kd' elements (one per variable in d())."
        exit 198
    }

    if `tau' <= 0 | `tau' >= 1 {
        di as error "tau() must lie in (0,1)."
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
    if "`distance'" == "mahalanobis" & `kx' < 2 {
        di as error "distance(mahalanobis) requires at least 2 covariates."
        exit 198
    }

    if "`method'" == "" local method "exact"
    if !inlist("`method'", "exact", "asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }

    if `s' == 0 {
        local s = floor(`n'^((`kx' + 1) / (`kx' + 2)))
    }
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

    local beta_mata ""
    foreach v of local beta {
        local beta_mata `beta_mata' `v',
    }
    local beta_mata = substr("`beta_mata'", 1, length("`beta_mata'") - 1)

    tempname b V beta_mat
    mata {
        _pl_X      = st_data(., "`xvars'", "`touse'")
        _pl_y      = st_data(., "`outcome'", "`touse'")
        _pl_D      = st_data(., "`dvars'", "`touse'")
        _pl_q      = st_data(., "`qvar'", "`touse'")
        _pl_md     = st_data(., "`mdvars'", "`touse'")
        _pl_x0     = (`x0_mata')
        _pl_beta   = (`beta_mata')'
        _pl_ranks  = rdistnn_ranks(_pl_X, _pl_x0, "`distance'")
        _pl_w      = rdistnn_weights_from_ranks(_pl_ranks, `k', `s', "`method'")
        _pl_resid  = `tau' :- (_pl_y :<= (_pl_D * _pl_beta :+ _pl_q))
        _pl_g      = (_pl_D :- _pl_md) :* (_pl_resid * J(1, cols(_pl_D), 1))
        _pl_moment = cross(_pl_w, _pl_g)
        _pl_V      = rdistnn_jk_cov_from_ranks(_pl_ranks, _pl_g, `k', `s', "`method'")
        _pl_neff   = 1 / sum(_pl_w :^2)
        st_matrix("`b'", _pl_moment)
        st_matrix("`V'", _pl_V)
        st_matrix("`beta_mat'", _pl_beta')
        st_numscalar("r(neff)", _pl_neff)
    }

    matrix rownames `b' = moment
    matrix colnames `b' = `dvars'
    matrix rownames `V' = `dvars'
    matrix colnames `V' = `dvars'
    matrix colnames `beta_mat' = `dvars'

    ereturn post `b' `V', esample(`touse')
    ereturn matrix beta = `beta_mat'
    ereturn scalar N      = `n'
    ereturn scalar neff   = scalar(r(neff))
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar kx     = `kx'
    ereturn scalar kd     = `kd'
    ereturn scalar tau    = `tau'
    ereturn local  x0       "`x0'"
    ereturn local  xvars    "`xvars'"
    ereturn local  dvars    "`dvars'"
    ereturn local  mdvars   "`mdvars'"
    ereturn local  qvar     "`qvar'"
    ereturn local  yvar     "`outcome'"
    ereturn local  distance "`distance'"
    ereturn local  method   "`method'"
    ereturn local  cmd      "plqr_dnn"
    ereturn local  cmdline  "plqr_dnn `0'"

    _plqr_dnn_Display, level(`level')
end

program define _plqr_dnn_Display
    syntax [, Level(cilevel)]
    di ""
    di as text "DNN PLQR Orthogonal Moments"
    di as text "Evaluation point x0: " as result "`e(x0)'"
    di as text "Covariates: "          as result "`e(xvars)'"
    di as text "Regressors: "          as result "`e(dvars)'"
    di as text "Distance: " as result "`e(distance)'" ///
       as text "  k = " as result e(k) ///
       as text "  s = " as result e(s) ///
       as text "  neff = " as result %6.2f e(neff)
    di as text "tau = " as result %6.3f e(tau)
    di ""
    ereturn display, level(`level')
end
