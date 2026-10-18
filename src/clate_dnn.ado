*! version 0.1.0  2026-10-17  Jakob Juergens
*! Experimental conditional LATE with supplied binary-IV nuisance predictions

program define clate_dnn, eclass
    version 17
    if replay() {
        if "`e(cmd)'" != "clate_dnn" error 301
        _clate_dnn_Display `0'
        exit
    }
    syntax varname(numeric) [if] [in], [ ///
        X(varlist numeric) XCOLS(varlist numeric) ///
        X0(numlist) X0MATRIX(name) ///
        DCOL(varname numeric) ACOL(varname numeric) ECOL(varname numeric) ///
        MY0COL(varname numeric) MY1COL(varname numeric) ///
        MD0COL(varname numeric) MD1COL(varname numeric) ///
        PSEUDOOutcome DENOMinator(varname numeric) ONESIDED ///
        CLIP(real 1e-8) DENOMINATORTOL(real 1.490116119384766e-8) ///
        TRAINING(string) K(integer 1) S(integer 0) ///
        TUNINGMode(string) LAMBDA(real .4) LAMBDAGrid(numlist) ///
        SGRID(numlist) SDEFAULT(integer 0) ///
        DISTance(string) METHod(string) Level(cilevel) noSE]

    local original_cmdline `"clate_dnn `0'"'
    local input `varlist'
    if "`x'"!="" & "`xcols'"!="" {
        di as error "Specify either x() or xcols(), not both."
        exit 198
    }
    local xvars `x' `xcols'
    if "`xvars'"=="" {
        di as error "xcols() is required (alias: x())."
        exit 198
    }
    if missing(`clip') | `clip'<=0 | `clip'>=.5 {
        di as error "clip() must be finite and strictly between zero and 0.5."
        exit 198
    }
    if missing(`denominatortol') | `denominatortol'<0 {
        di as error "denominatortol() must be finite and nonnegative."
        exit 198
    }
    if "`training'"=="" local training "supplied"
    if !inlist("`training'","supplied","oracle","full_sample","oob","cross_fit") {
        di as error "training() must be supplied, oracle, full_sample, oob, or cross_fit."
        exit 198
    }
    local route "augmented_iv"
    local nuisance_vars `dcol' `acol' `ecol' `my0col' `my1col' `md0col' `md1col'
    if "`pseudooutcome'"!="" {
        local route "paired_signals"
        if "`denominator'"=="" | "`nuisance_vars'`onesided'"!="" {
            di as error "pseudooutcome requires denominator() and cannot be combined with nuisance options or onesided."
            exit 198
        }
    }
    else {
        if "`denominator'"!="" {
            di as error "denominator() requires pseudooutcome."
            exit 198
        }
        foreach arg in dcol acol ecol my0col my1col md1col {
            if "``arg''"=="" {
                di as error "`arg'() is required for augmented IV signals."
                exit 198
            }
        }
        if "`md0col'"=="" & "`onesided'"=="" {
            di as error "md0col() is required unless onesided is specified."
            exit 198
        }
    }

    marksample touse, novarlist
    if "`route'"=="augmented_iv" {
        foreach prob in ecol md0col md1col {
            if "``prob''"!="" {
                quietly count if `touse' & (missing(``prob'') | ``prob''<0 | ``prob''>1)
                if r(N)>0 {
                    di as error "`prob'() must contain finite predictions in [0,1] before clipping."
                    exit 198
                }
            }
        }
        markout `touse' `input' `xvars' `nuisance_vars'
        foreach binary in dcol acol {
            quietly count if `touse' & !inlist(``binary'',0,1)
            if r(N)>0 {
                di as error "`binary'() must be coded 0/1."
                exit 198
            }
        }
        if "`onesided'"!="" {
            quietly count if `touse' & `dcol'>`acol'
            if r(N)>0 {
                di as error "onesided requires D <= A on the estimation sample."
                exit 198
            }
            if "`md0col'"!="" {
                quietly count if `touse' & `md0col'!=0
                if r(N)>0 {
                    di as error "md0col() must be zero when onesided is specified."
                    exit 198
                }
            }
        }
    }
    else markout `touse' `input' `denominator' `xvars'
    quietly count if `touse'
    if r(N)==0 error 2000

    tempvar phiy phid eclip
    if "`route'"=="augmented_iv" {
        local md0value `md0col'
        if "`md0value'"=="" local md0value "0"
        quietly gen double `eclip' = min(max(`ecol',`clip'),1-`clip') if `touse'
        quietly gen double `phiy' = `my1col'-`my0col' + ///
            `acol'*(`input'-`my1col')/`eclip' - ///
            (1-`acol')*(`input'-`my0col')/(1-`eclip') if `touse'
        quietly gen double `phid' = `md1col'-`md0value' + ///
            `acol'*(`dcol'-`md1col')/`eclip' - ///
            (1-`acol')*(`dcol'-`md0value')/(1-`eclip') if `touse'
    }
    else {
        quietly gen double `phiy' = `input' if `touse'
        quietly gen double `phid' = `denominator' if `touse'
    }
    quietly count if `touse' & (missing(`phiy') | missing(`phid'))
    if r(N)>0 {
        di as error "Inputs produce non-finite augmented signals; rescale outcome inputs."
        exit 498
    }

    // Use the existing CATE localization/tuning contract, then replace its
    // scalar-signal calculation with the paired ratio calculation below.
    local opts "xcols(`xvars') k(`k') s(`s') lambda(`lambda') sdefault(`sdefault') level(`level')"
    foreach arg in x0 x0matrix tuningmode lambdagrid sgrid distance method {
        if "``arg''"!="" local opts "`opts' `arg'(``arg'')"
    }
    // The shared jackknife short-circuits zero signals, so this geometry/tuning
    // delegation performs no unused scalar deletion/covariance calculation.
    // Keep SE feasibility rules active when paired inference is requested.
    tempvar localization_signal
    quietly gen byte `localization_signal' = 0 if `touse'
    quietly cate_dnn `localization_signal' if `touse', pseudooutcome `opts' `se'
    local n = e(N)
    local m = e(m)
    local resolved_s = e(s)
    local resolved_method "`e(method)'"
    local resolved_distance "`e(distance)'"
    local isgrid = ("`e(tuning_mode)'"=="grid")
    local active = ("`se'"=="")
    tempname targets diag components theta vcov componentV status semat ///
        meta gridtheta gridse gridy gridd gridstatus
    matrix `targets' = e(targets)
    matrix `diag' = e(diagnostics)
    if `isgrid' matrix `meta' = e(tuning_candidates)
    else matrix `meta' = (1,1,1,1,0,1,.,`resolved_s')
    local ng = rowsof(`meta')
    local default = 1
    if `isgrid' local default = e(tuning_default_candidate)

    findfile "_rdistnn_clate_mata.ado"
    quietly run "`r(fn)'"
    mata {
        _cl_X = st_data(., "`xvars'", "`touse'")
        _cl_X0 = st_matrix("`targets'")
        _cl_ranks = rdistnn_ranks_matrix(_cl_X,_cl_X0,"`resolved_distance'")
        _cl_signals = st_data(., "`phiy' `phid'", "`touse'")
        _cl_meta = st_matrix("`meta'")
        _cl_grid = _cl_gridse = _cl_gridy = _cl_gridd = J(`ng',`m',.)
        _cl_gridstatus = J(`ng'*`m',7,.)
        _cl_c = _cl_V = _cl_CV = _cl_status = J(0,0,.)
        _cl_b = J(1,0,.)
        for (_cl_g=1; _cl_g<=`ng'; _cl_g++) {
            if (_cl_meta[_cl_g,4]==0) continue
            _cl_prior=0
            for (_cl_h=1; _cl_h<_cl_g; _cl_h++) {
                if (_cl_meta[_cl_h,4] & _cl_meta[_cl_h,2]==_cl_meta[_cl_g,2]) {
                    _cl_prior=_cl_h
                    break
                }
            }
            if (_cl_prior) {
                _cl_grid[_cl_g,.]=_cl_grid[_cl_prior,.]
                _cl_gridse[_cl_g,.]=_cl_gridse[_cl_prior,.]
                _cl_gridy[_cl_g,.]=_cl_gridy[_cl_prior,.]
                _cl_gridd[_cl_g,.]=_cl_gridd[_cl_prior,.]
                _cl_gridstatus[((_cl_g-1)*`m'+1)..(_cl_g*`m'),.] =
                    (J(`m',1,_cl_g),_cl_gridstatus[((_cl_prior-1)*`m'+1)..(_cl_prior*`m'),2..7])
                continue
            }
            rdistnn_clate_fit(_cl_ranks,_cl_signals,`k',_cl_meta[_cl_g,8],
                "`resolved_method'",`active',`denominatortol',
                _cl_c,_cl_b,_cl_V,_cl_CV,_cl_status)
            _cl_grid[_cl_g,.] = _cl_b
            _cl_gridse[_cl_g,.] = sqrt(diagonal(_cl_V))'
            _cl_gridy[_cl_g,.] = _cl_c[,1]'
            _cl_gridd[_cl_g,.] = _cl_c[,2]'
            _cl_gridstatus[((_cl_g-1)*`m'+1)..(_cl_g*`m'),.] =
                (J(`m',1,_cl_g),(1::`m'),_cl_status)
            if (_cl_meta[_cl_g,2]==_cl_meta[`default',2]) {
                st_matrix("`components'",_cl_c)
                st_matrix("`theta'",_cl_b)
                st_matrix("`vcov'",_cl_V)
                st_matrix("`componentV'",_cl_CV)
                st_matrix("`status'",_cl_status)
            }
        }
        _cl_b = st_matrix("`theta'")
        _cl_W = rdistnn_weights_from_rmat(_cl_ranks,`k',`resolved_s',"`resolved_method'")
        _cl_A = _cl_e = J(0,1,.)
        if ("`route'"=="augmented_iv") {
            _cl_A = st_data(.,"`acol'","`touse'")
            _cl_e = st_data(.,"`eclip'","`touse'")
        }
        st_matrix("`diag'",rdistnn_clate_diagnostics(st_matrix("`diag'"),
            _cl_W,_cl_signals,_cl_b,_cl_A,_cl_e))
        st_matrix("`semat'",sqrt(diagonal(st_matrix("`vcov'")))')
        st_matrix("`gridtheta'",_cl_grid)
        st_matrix("`gridse'",_cl_gridse)
        st_matrix("`gridy'",_cl_gridy)
        st_matrix("`gridd'",_cl_gridd)
        st_matrix("`gridstatus'",_cl_gridstatus)
    }

    // Preserve the shared localization results across a fresh, honest post.
    // Stata forbids missing coefficients: undefined targets are omitted from
    // e(b), while all target identities remain in the canonical full tables.
    local saved_scalars : e(scalars)
    local saved_macros : e(macros)
    local saved_matrices : e(matrices)
    foreach name of local saved_scalars {
        local ss_`name' = e(`name')
    }
    foreach name of local saved_macros {
        local sl_`name' `"`e(`name')'"'
    }
    foreach name of local saved_matrices {
        if inlist("`name'","b","V") continue
        tempname sm_`name'
        matrix `sm_`name'' = e(`name')
    }
    local cnames ""
    local validnames ""
    local keep ""
    local compnames ""
    forvalues j=1/`m' {
        local coef "theta`j'"
        if `m'==1 local coef "theta"
        local cnames `cnames' `coef'
        local compnames `compnames' y`j' d`j'
        if el(`theta',1,`j')<. {
            local keep `keep' `j'
            local validnames `validnames' `coef'
        }
    }
    if "`keep'"!="" {
        tempname postb postV
        mata: _cl_keep=strtoreal(tokens("`keep'")); st_matrix("`postb'",st_matrix("`theta'")[,_cl_keep]); st_matrix("`postV'",st_matrix("`vcov'")[_cl_keep,_cl_keep]); st_local("postcov",strofreal(!hasmissing(st_matrix("`postV'"))))
        matrix colnames `postb' = `validnames'
        matrix colnames `postV' = `validnames'
        matrix rownames `postV' = `validnames'
        if `active' & `postcov' ereturn post `postb' `postV', esample(`touse')
        else ereturn post `postb', esample(`touse')
    }
    else ereturn post, esample(`touse')
    foreach name of local saved_scalars {
        ereturn scalar `name' = `ss_`name''
    }
    foreach name of local saved_macros {
        if "`name'"=="properties" continue
        ereturn local `name' `"`sl_`name''"'
    }
    foreach name of local saved_matrices {
        if inlist("`name'","b","V") continue
        ereturn matrix `name' = `sm_`name'', copy
    }
    matrix colnames `theta' = `cnames'
    matrix colnames `semat' = `cnames'
    matrix colnames `vcov' = `cnames'
    matrix rownames `vcov' = `cnames'
    matrix colnames `componentV' = `compnames'
    matrix rownames `componentV' = `compnames'
    matrix rownames `components' = `: rownames `targets''
    matrix colnames `components' = reduced_form first_stage
    matrix rownames `status' = `: rownames `targets''
    matrix colnames `status' = ratio_status negative_denominator ///
        near_zero_denominator covariance_status wald_available
    matrix rownames `diag' = `: rownames `targets''
    matrix colnames `diag' = `: colnames `sm_diagnostics''
    ereturn matrix estimates = `theta', copy
    ereturn matrix se_targets = `semat', copy
    ereturn matrix theta_vcov = `vcov', copy
    ereturn matrix component_vcov = `componentV', copy
    ereturn matrix components = `components', copy
    ereturn matrix status = `status', copy
    ereturn matrix diagnostics = `diag', copy
    ereturn matrix target_diagnostics = `diag', copy
    if `isgrid' {
        foreach mat in gridtheta gridse gridy gridd {
            matrix colnames ``mat'' = `cnames'
            matrix rownames ``mat'' = `: rownames `meta''
        }
        matrix colnames `gridstatus' = candidate_id target_id ratio_status ///
            negative_denominator near_zero_denominator covariance_status wald_available
        ereturn matrix tuning_estimates = `gridtheta', copy
        ereturn matrix tuning_se = `gridse', copy
        ereturn matrix tuning_reduced_form = `gridy', copy
        ereturn matrix tuning_first_stage = `gridd', copy
        ereturn matrix tuning_target_status = `gridstatus', copy
    }
    ereturn scalar clip = `clip'
    ereturn scalar denominator_tol = `denominatortol'
    ereturn scalar one_sided = ("`onesided'"!="")
    ereturn local cmd "clate_dnn"
    ereturn local cmdline `"`original_cmdline'"'
    ereturn local score_method "`route'"
    ereturn local se_method = cond(`active',"paired_rank_shift","none")
    ereturn local covariance = cond(`active',"paired_frozen_score","none")
    ereturn local outcome "`input'"
    ereturn local signal "`input'"
    ereturn local denominator "`denominator'"
    ereturn local d_col "`dcol'"
    ereturn local a_col "`acol'"
    ereturn local e_col "`ecol'"
    ereturn local treat "`acol'"
    ereturn local propensity "`ecol'"
    ereturn local propensity_version = cond("`route'"=="augmented_iv","clipped","not_used")
    ereturn local my0_col "`my0col'"
    ereturn local my1_col "`my1col'"
    ereturn local md0_col "`md0col'"
    ereturn local md1_col "`md1col'"
    ereturn local training "`training'"
    ereturn local overlap_role "instrument"
    ereturn local ratio_status_codes "0=ok 1=zero_denominator 2=numerical_failure"
    ereturn local covariance_status_codes "0=ok 1=not_computed 2=zero_score_variance 3=numerical_failure"
    local signature_vars `xvars' `input' `nuisance_vars' `denominator'
    local signature_vars : list uniq signature_vars
    quietly _rdistnn_estimation_state `signature_vars'
    tempname model_signals model_ranks model_distances
    mata: st_matrix("`model_signals'",_cl_signals); st_matrix("`model_ranks'",_cl_ranks); st_matrix("`model_distances'",rdistnn_distance_matrix(_cl_X,_cl_X0,"`resolved_distance'"))
    quietly _rdistnn_model_store `signature_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances')
    _clate_dnn_Display, level(`level')
end

program define _clate_dnn_Display, rclass
    version 17
    _rdistnn_report `0'
    return add
end
