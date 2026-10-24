*! version 0.1.0  2026-07-12  Jakob Juergens
*! Experimental PLQR localized-moment evaluation at supplied theta

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
        MD(varlist numeric)              ///  supplied m_D(X), same order as d()
        G(varname numeric)               ///  supplied g_tau(X)
        THETAEVAL(numlist)                ///  supplied theta vector
        [                                ///
        tau(real 0.5)                    ///  quantile index
        k(integer 1)                     ///  # selected within each s-subset
        s(integer 0)                     ///  subset/kernel order; larger = tighter
        Distance(string)                 ///  localization geometry
        SCALERule(string) SCALEVec(numlist) SCALEEPS(real 1e-12) ///
        Method(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///
        ]

    local gvar `g'
    local theta `thetaeval'

    marksample touse
    markout `touse' `varlist' `x' `d' `gvar' `md'

    local outcome `varlist'
    local xvars   `x'
    local dvars   `d'
    local mdvars  `md'
    local gvar    `gvar'
    local kx : word count `xvars'
    local kd : word count `dvars'
    local kmd : word count `mdvars'
    local ntheta : word count `theta'

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

    if `ntheta' != `kd' {
        di as error "thetaeval() must have `kd' elements (one per variable in d())."
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
    if !inlist("`distance'", "l2", "mahalanobis", "l1", "l1_scaled", "l2_scaled") {
        di as error "distance() must be l1, l2, l1_scaled, l2_scaled, or mahalanobis."
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
        local s = floor(.4 * `n'^((`kx' + 1) / (`kx' + 2)))
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

    tempname b V theta_mat moment_mat moment_vcov_mat
    tempname targets model_groups model_signals model_ranks model_distances model_geometry diagnostics
    mata: rdistnn_plqr_fit_ado()

    matrix rownames `b' = moment
    matrix colnames `b' = `dvars'
    matrix rownames `V' = `dvars'
    matrix colnames `V' = `dvars'
    matrix colnames `theta_mat' = `dvars'
    matrix `moment_mat' = `b'
    matrix `moment_vcov_mat' = `V'

    ereturn post `b' `V', esample(`touse')
    ereturn matrix moment = `moment_mat', copy
    ereturn matrix moment_vcov = `moment_vcov_mat', copy
    ereturn matrix theta_eval = `theta_mat', copy
    ereturn scalar N      = `n'
    ereturn scalar neff   = scalar(r(neff))
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar kx     = `kx'
    ereturn scalar kd     = `kd'
    ereturn scalar m      = 1
    ereturn scalar p      = `kx'
    ereturn scalar tau    = `tau'
    if "`scalerule'"=="" local scalerule "sd"
    ereturn local scalerule "`scalerule'"
    ereturn scalar scaleeps=`scaleeps'
    ereturn local  x0       "`x0'"
    ereturn local  xvars    "`xvars'"
    ereturn local  dvars    "`dvars'"
    ereturn local  mdvars   "`mdvars'"
    ereturn local  gvar     "`gvar'"
    ereturn local  theta_source "supplied"
    ereturn local  inference_target "moment"
    ereturn local  root_solver "none"
    ereturn local  yvar     "`outcome'"
    ereturn local  distance "`distance'"
    ereturn local  method   "`method'"
    ereturn local  se_method "rank_shift"
    ereturn local  jackknife "fast"
    ereturn local  tie_rule "exact_distance_block_averaging_v1"
    ereturn local  cmd      "plqr_dnn"
    ereturn local  cmdline  "plqr_dnn `0'"

    tempname ses
    matrix colnames `targets'=`xvars'
    ereturn matrix targets=`targets', copy
    local snapshot_vars `xvars' `outcome' `dvars' `gvar' `mdvars'
    local snapshot_vars : list uniq snapshot_vars
    local geometry_option ""
    if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        local geometry_option "geometry(`model_geometry')"
    }
    mata: st_matrix("`ses'",sqrt(diagonal(st_matrix("e(V)")))')
    matrix colnames `ses'=`dvars'
    ereturn matrix estimates=`moment_mat', copy
    ereturn matrix se_targets=`ses', copy
    ereturn matrix estimate_vcov=`moment_vcov_mat', copy
    tempname bytarget
    matrix `bytarget'=`moment_mat''
    ereturn matrix moment_by_target=`bytarget', copy
    foreach file in _rdistnn_clate_mata.ado _rdistnn_model_mata.ado {
        findfile "`file'"
        quietly run "`r(fn)'"
    }
    matrix colnames `diagnostics'=ess max_weight top1_weight_share top5_weight_share top10_weight_share n_nonzero_weight distance_available radius_rank_k radius_rank_s weighted_mean_distance weighted_median_distance
    ereturn matrix diagnostics=`diagnostics', copy
    ereturn matrix target_diagnostics=`diagnostics', copy
    quietly _rdistnn_model_store `snapshot_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances') groups(`model_groups') `geometry_option'
    quietly signestimationsample `snapshot_vars'
    quietly _rdistnn_sample_signature `snapshot_vars' if e(sample)
    ereturn local ordered_signature "`r(signature)'"
    ereturn local signature_vars "`snapshot_vars'"
    _plqr_dnn_Display, level(`level')
end

program define _plqr_dnn_Display, rclass
    syntax [, Level(cilevel)]
    _rdistnn_report, level(`level')
    return add
end
