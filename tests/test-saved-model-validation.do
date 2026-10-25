// Corrupt current snapshots must fail before reuse or observation reconstruction.
version 17
do "tests/test_framework.do"
do "tests/helpers/e-state.do"
local smv_distance "$rdistnn_distance_engine"
local smv_jackknife "$rdistnn_jackknife_engine"
global rdistnn_distance_engine mata
global rdistnn_jackknife_engine mata
clear
set obs 8
gen double smv_x=_n-4.5
gen double smv_z=0
gen double smv_y=_n^2
gen double smv_den=1+_n/100
gen double smv_d=2+sin(_n)
gen double smv_d2=-3+cos(_n)
gen double smv_zero=0
matrix SMV_TARGETS=(0,0\2.5,0)
quietly np_dnn smv_y, xcols(smv_x smv_z) x0matrix(SMV_TARGETS) s(3)
estimates store SMV_PARENT
capture program drop _smv_corrupt
program define _smv_corrupt, eclass
    syntax, CASE(string)
    tempname value
    if inlist("`case'","negative_distance","missing_distance") {
        matrix `value'=e(model_distances)
        matrix `value'[1,1]=cond("`case'"=="negative_distance",-1,.)
        ereturn matrix model_distances=`value', copy
    }
    if inlist("`case'","rank_order","rank_tie","duplicate_rank") {
        matrix `value'=e(model_ranks)
        local a=1
        local b=2
        if "`case'"=="rank_tie" {
            local a=4
            local b=5
        }
        local rank=el(`value',`a',1)
        matrix `value'[`a',1]=el(`value',`b',1)
        if "`case'"!="duplicate_rank" matrix `value'[`b',1]=`rank'
        ereturn matrix model_ranks=`value', copy
    }
    if inlist("`case'","fractional_id","duplicate_id","missing_id","short_ids") {
        matrix `value'=e(model_row_ids)
        if "`case'"=="fractional_id" matrix `value'[1,1]=-1.5
        if "`case'"=="duplicate_id" matrix `value'[2,1]=el(`value',1,1)
        if "`case'"=="missing_id" matrix `value'[1,1]=.
        if "`case'"=="short_ids" matrix `value'=`value'[1..7,1]
        ereturn matrix model_row_ids=`value', copy
    }
    if "`case'"=="input_rows" {
        matrix `value'=e(model_inputs)
        matrix `value'=`value'[1..7,1...]
        ereturn matrix model_inputs=`value', copy
    }
    if "`case'"=="input_names" ereturn local model_vars "smv_x smv_z"
    if "`case'"=="signal_columns" {
        matrix `value'=e(model_signals), e(model_signals)
        ereturn matrix model_signals=`value', copy
    }
    if "`case'"=="signal_missing" {
        matrix `value'=e(model_signals)
        matrix `value'[1,1]=.
        ereturn matrix model_signals=`value', copy
    }
    if "`case'"=="method" ereturn local method "invalid"
    if "`case'"=="distance" ereturn local distance "invalid"
    if "`case'"=="count" ereturn scalar N=9
    if "`case'"=="dimension" ereturn scalar p=3
    if "`case'"=="target_count" ereturn scalar m=3
    if "`case'"=="legacy" ereturn scalar model_version=2
    if "`case'"=="old_policy" ereturn local tie_rule "distance_then_covariates_with_duplicate_block_averaging"
    if "`case'"=="group_mismatch" {
        matrix `value'=e(model_groups)
        matrix `value'[1,1]=99
        ereturn matrix model_groups=`value', copy
    }
    if "`case'"=="group_missing" {
        matrix `value'=(.)
        ereturn matrix model_groups=`value', copy
    }
    if "`case'"=="geometry" {
        matrix `value'=(1,0)
        ereturn matrix model_geometry=`value', copy
    }
    if inlist("`case'","target_identity","target_dimensions") {
        matrix `value'=e(targets)
        if "`case'"=="target_identity" matrix `value'[1,1]=987654321
        else matrix `value'=`value'[1,1...]
        ereturn matrix targets=`value', copy
    }
    if inlist("`case'","estimates_dimensions","ses_dimensions","vcov_dimensions","diagnostics_dimensions") {
        local field "estimates"
        if "`case'"=="ses_dimensions" local field "se_targets"
        if "`case'"=="vcov_dimensions" local field "estimate_vcov"
        if "`case'"=="diagnostics_dimensions" local field "diagnostics"
        matrix `value'=(1)
        ereturn matrix `field'=`value', copy
    }
    if "`case'"=="missing_estimate" {
        matrix `value'=e(estimates)
        matrix `value'[1,1]=.
        ereturn matrix estimates=`value', copy
    }
    if inlist("`case'","status_dimensions","components_dimensions","componentcov_dimensions") {
        local field "status"
        if "`case'"=="components_dimensions" local field "components"
        if "`case'"=="componentcov_dimensions" local field "component_vcov"
        matrix `value'=(1)
        ereturn matrix `field'=`value', copy
    }
    if inlist("`case'","components_missing","components_zero_denominator") {
        matrix `value'=e(components)
        if "`case'"=="components_missing" matrix `value'[1,1]=.
        else matrix `value'[1,2]=0
        ereturn matrix components=`value', copy
    }
    if "`case'"=="component_sign_flag" {
        matrix `value'=e(status)
        matrix `value'[1,2]=1-el(`value',1,2)
        ereturn matrix status=`value', copy
    }
    if substr("`case'",1,13)=="componentcov_" {
        matrix `value'=e(component_vcov)
        if "`case'"=="componentcov_negative" matrix `value'[1,1]=-1
        if "`case'"=="componentcov_asymmetric" matrix `value'[1,2]=el(`value',1,2)+1
        if "`case'"=="componentcov_partial" matrix `value'[1,1]=.
        if "`case'"=="componentcov_cross_missing" {
            matrix `value'[1,3]=.
            matrix `value'[3,1]=.
        }
        if "`case'"=="componentcov_false_success" {
            forvalues j=1/`=colsof(`value')' {
                forvalues h=1/2 {
                    matrix `value'[`h',`j']=.
                    matrix `value'[`j',`h']=.
                }
            }
        }
        if "`case'"=="componentcov_nose_finite" matrix `value'[1,1]=0
        ereturn matrix component_vcov=`value', copy
    }
    if "`case'"=="negative_se" {
        matrix `value'=e(se_targets)
        matrix `value'[1,1]=-2
        ereturn matrix se_targets=`value', copy
    }
    if inlist("`case'","missing_se","inconsistent_se") {
        matrix `value'=e(se_targets)
        matrix `value'[1,1]=cond("`case'"=="missing_se",.,el(`value',1,1)+1)
        ereturn matrix se_targets=`value', copy
    }
    if inlist("`case'","negative_variance","asymmetric_vcov","missing_variance") {
        matrix `value'=e(estimate_vcov)
        if "`case'"=="negative_variance" matrix `value'[1,1]=-1
        if "`case'"=="missing_variance" matrix `value'[1,1]=.
        if "`case'"=="asymmetric_vcov" matrix `value'[1,2]=el(`value',1,2)+1
        ereturn matrix estimate_vcov=`value', copy
    }
    if "`case'"=="se_method" ereturn local se_method "invalid"
    if "`case'"=="jackknife_method" ereturn local jackknife "none"
    if inlist("`case'","status_values","status_missing","status_availability") {
        matrix `value'=e(status)
        if "`case'"=="status_values" matrix `value'[1,5]=2
        if "`case'"=="status_missing" matrix `value'[1,1]=.
        if "`case'"=="status_availability" matrix `value'[1,1]=1
        ereturn matrix status=`value', copy
    }
    if inlist("`case'","status_ok_missing","status_ok_zero","status_zero_nonzero","status_failed_finite") {
        tempname ses covariance
        matrix `value'=e(status)
        matrix `ses'=e(se_targets)
        matrix `covariance'=e(estimate_vcov)
        matrix `value'[1,5]=0
        if inlist("`case'","status_ok_missing","status_ok_zero") {
            local v=cond("`case'"=="status_ok_missing",.,0)
            matrix `ses'[1,1]=`v'
            forvalues j=1/`=colsof(`covariance')' {
                matrix `covariance'[1,`j']=`v'
                matrix `covariance'[`j',1]=`v'
            }
        }
        if "`case'"=="status_zero_nonzero" matrix `value'[1,4]=2
        if "`case'"=="status_failed_finite" matrix `value'[1,4]=3
        ereturn matrix status=`value', copy
        ereturn matrix se_targets=`ses', copy
        ereturn matrix estimate_vcov=`covariance', copy
    }
    if "`case'"=="bad_arm_link" ereturn local t_col "not_retained"
    if "`case'"=="bad_propensity_link" ereturn local e_col "not_retained"
    if "`case'"=="missing_arm_link" ereturn local t_col ""
    if "`case'"=="missing_propensity_link" ereturn local e_col ""
    if "`case'"=="negative_clip" ereturn scalar clip=-.2
    if "`case'"=="missing_clip" ereturn scalar clip=.
    if "`case'"=="large_clip" ereturn scalar clip=.9
    if "`case'"=="negative_tolerance" ereturn scalar denominator_tol=-1
    if "`case'"=="missing_tolerance" ereturn scalar denominator_tol=.
    if "`case'"=="signature_vars_missing" ereturn local signature_vars ""
    if "`case'"=="signature_vars_mismatch" ereturn local signature_vars "wrong retained signature variables"
    if "`case'"=="ordered_signature_missing" ereturn local ordered_signature ""
    if substr("`case'",1,5)=="grid_" {
        local field "tuning_candidates"
        if "`case'"=="grid_estimates_shape" local field "tuning_estimates"
        if "`case'"=="grid_se_shape" local field "tuning_se"
        if "`case'"=="grid_diag_shape" local field "tuning_diagnostics"
        if "`case'"=="grid_status_shape" local field "tuning_target_status"
        if "`case'"=="grid_y_shape" local field "tuning_reduced_form"
        if "`case'"=="grid_d_shape" local field "tuning_first_stage"
        if strpos("`case'","_shape") {
            matrix `value'=(1)
            ereturn matrix `field'=`value', copy
            exit
        }
        if "`case'"=="grid_alias" {
            matrix `value'=e(tuning_estimates)
            matrix `value'[3,1]=el(`value',3,1)+1
            ereturn matrix tuning_estimates=`value', copy
            exit
        }
        if "`case'"=="grid_availability" {
            matrix `value'=e(tuning_estimates)
            matrix `value'[4,1]=1
            ereturn matrix tuning_estimates=`value', copy
            exit
        }
        if "`case'"=="grid_target_names" {
            matrix `value'=e(tuning_estimates)
            local names ""
            forvalues j=1/`=e(m)' {
                local names "`names' wrong`j'"
            }
            matrix colnames `value'=`names'
            ereturn matrix tuning_estimates=`value', copy
            exit
        }
        if "`case'"=="grid_count" {
            ereturn scalar tuning_n_unique=e(tuning_n_unique)+1
            exit
        }
        matrix `value'=e(tuning_candidates)
        if "`case'"=="grid_default" matrix `value'[1,3]=1
        if "`case'"=="grid_id" matrix `value'[1,1]=77
        if "`case'"=="grid_computed" matrix `value'[1,2]=77
        if "`case'"=="grid_feasible" matrix `value'[1,4]=0
        ereturn matrix tuning_candidates=`value', copy
    }
end

// Both retained-state fixtures exercise the same public routes and e() guard.
capture program drop _smv_check_routes
program define _smv_check_routes
    syntax, ESTIMATES(name) SNAPSHOT(name) CASE(string) ROUTES(string) [FAMILY(name) QUERY CALLERTEST(string)]
    if "`family'"=="" local family "np"
    local label "`case'"
    local rejected "corrupt snapshot is rejected"
    local preserved "failure preserves complete e() and sample"
    local caller "_N==8 & smv_y[1]==1"
    local graph "smv_invalid"
    if "`query'"!="" {
        local label "`family'/`case'"
        local rejected "reject invalid retained state with refit error"
        local preserved "preserve current e() and query data"
        local caller `"_N==3 & caller[1]=="preserved" & smvc_x[1]==777"'
        local graph "smvc_invalid"
    }
    if `"`callertest'"'!="" local caller `"`callertest'"'
    foreach route of local routes {
        quietly estimates restore `estimates'
        if "`route'"=="cached" capture quietly rdistnn_predict, x0matrix(SMV_TARGETS) post
        if "`route'"=="new" capture quietly rdistnn_predict, x0(.25 0) post
        if "`route'"=="replay" capture quietly `family'_dnn
        if inlist("`route'","estimates","results") capture quietly rdistnn_results
        if "`route'"=="vcov" capture quietly rdistnn_results, type(vcov)
        if "`route'"=="components" capture quietly rdistnn_results, type(components)
        if "`route'"=="resultdiag" capture quietly rdistnn_results, type(diagnostics)
        if "`route'"=="diagnostics" capture quietly `family'_dnn_diagnostics
        if "`route'"=="weights" capture quietly rdistnn_results, type(weights)
        if "`route'"=="kernel" capture quietly `family'_dnn_plot, type(weights) view(kernel) name(`graph') replace nodraw
        if "`route'"=="profile" capture quietly `family'_dnn_plot, type(weights) name(`graph') replace nodraw
        if "`route'"=="estimateplot" capture quietly `family'_dnn_plot, type(estimates) name(`graph') replace nodraw
        if "`route'"=="diagplot" capture quietly `family'_dnn_plot, type(diagnostics) name(`graph') replace nodraw
        if "`route'"=="plot" capture quietly `family'_dnn_plot, name(`graph') replace nodraw
        if "`route'"=="predict" capture quietly rdistnn_predict, x0(.25 0) sefit post
        if "`route'"=="candidate" capture quietly rdistnn_candidate, candidate(1) post
        local rc=_rc
        test_assert "`label'/`route': `rejected'" `rc'==498
        _rdistnn_test_e_compare, snapshot(`snapshot')
        local same=r(equal)
        test_assert "`label'/`route': `preserved'" `same' & `caller'
    }
end

test_section "saved snapshot corruption: reject consistently and preserve caller state"
foreach case in negative_distance missing_distance rank_order rank_tie duplicate_rank ///
    fractional_id duplicate_id missing_id short_ids input_rows input_names ///
    signature_vars_missing signature_vars_mismatch ordered_signature_missing ///
    signal_columns signal_missing method distance count dimension target_count legacy old_policy group_mismatch group_missing {
    quietly estimates restore SMV_PARENT
    _smv_corrupt, case(`case')
    estimates store SMV_CORRUPT_ACTIVE
    quietly estimates restore SMV_CORRUPT_ACTIVE
    _rdistnn_test_e_snapshot, name(SMV_CORRUPT)
    _smv_check_routes, estimates(SMV_CORRUPT_ACTIVE) snapshot(SMV_CORRUPT) case(`case') ///
        routes("cached new weights kernel profile estimates vcov resultdiag diagnostics estimateplot diagplot")
    estimates drop SMV_CORRUPT_ACTIVE
    _rdistnn_test_e_snapshot_drop, name(SMV_CORRUPT)
}

// Malformed frozen geometry is invalid even when the requested targets are cached.
quietly np_dnn smv_y, xcols(smv_x smv_z) x0matrix(SMV_TARGETS) s(3) distance(l2_scaled)
_smv_corrupt, case(geometry)
capture quietly rdistnn_predict, x0matrix(SMV_TARGETS)
local rc=_rc
test_assert "cached targets reject a nonpositive frozen scale" `rc'==498
capture quietly rdistnn_results, type(weights)
local rc=_rc
test_assert "weight extraction rejects a nonpositive frozen scale" `rc'==498

test_section "valid cold snapshots: all estimator families retain fitted geometry and weights"
replace smv_x=smv_x[1] in 2
foreach family in np tdnn cate clate plqr {
    local covariates "xcols(smv_x smv_z) x0matrix(SMV_TARGETS)"
    if "`family'"=="np" quietly np_dnn smv_y, `covariates' s(3) distance(l2_scaled)
    if "`family'"=="tdnn" quietly np_tdnn smv_y, `covariates' s1(2) s2(3) distance(l2_scaled)
    if "`family'"=="cate" quietly cate_dnn smv_y, pseudooutcome `covariates' s(3) distance(l2_scaled)
    if "`family'"=="clate" quietly clate_dnn smv_y, pseudooutcome denominator(smv_den) `covariates' s(3) distance(l2_scaled)
    if "`family'"=="plqr" quietly plqr_dnn smv_y, x(smv_x smv_z) x0(0 0) d(smv_d smv_d2) md(smv_zero smv_zero) g(smv_zero) thetaeval(0 0) s(3) distance(l2_scaled)
    matrix SMV_FIT=e(estimates)
    matrix SMV_GEOMETRY=e(model_geometry)
    matrix SMV_CURRENT=e(model_targets)
    tempfile cold
    quietly estimates save "`cold'", replace
    preserve
    clear
    mata: mata clear
    quietly estimates use "`cold'"
    quietly rdistnn_predict, x0matrix(SMV_CURRENT) sefit
    matrix SMV_GOT=r(fit)
    mata: st_numscalar("smv_same",all(vec(st_matrix("SMV_GOT"))':==st_matrix("SMV_FIT")) & all(st_matrix("SMV_GEOMETRY"):==st_matrix("e(model_geometry)")))
    test_assert "`family': cold cached prediction keeps estimates and frozen geometry" scalar(smv_same) & _N==0
    quietly rdistnn_results, type(weights)
    local rows=rowsof(r(data))
    test_assert "`family': valid cold observation extraction keeps every ID" `rows'==e(N)*e(m)
    quietly rdistnn_predict, x0(.25 0) sefit post
    local columns=colsof(r(fit))
    test_assert "`family': valid cold posted prediction supports new targets" `columns'==1 & _N==0 & "`e(signature_vars)'"=="`e(model_vars)'" & "`e(ordered_signature)'"!=""
    foreach case in signature_vars_missing signature_vars_mismatch ordered_signature_missing {
        quietly estimates use "`cold'"
        _smv_corrupt, case(`case')
        capture quietly rdistnn_predict, x0(.25 0) post
        local rc=_rc
        test_assert "`family': missing or inconsistent `case' requests refitting" `rc'==498 & _N==0
    }
    quietly estimates use "`cold'"
    if inlist("`family'","cate","clate") {
        quietly rdistnn_results, type(influence)
        local rows=rowsof(r(data))
        test_assert "`family': valid cold influence retains every row" `rows'==e(N)*e(m)
        _smv_corrupt, case(negative_distance)
        capture quietly rdistnn_results, type(influence)
        local rc=_rc
        test_assert "`family': influence rejects corrupt saved distances" `rc'==498
    }
    if inlist("`family'","tdnn","plqr","clate") {
        quietly estimates use "`cold'"
        _smv_corrupt, case(signal_columns)
        capture quietly rdistnn_predict, x0matrix(SMV_CURRENT)
        local rc=_rc
        test_assert "`family': wrong signal component count is rejected" `rc'==498
    }
    restore
}
quietly estimates restore SMV_PARENT
estimates drop SMV_PARENT
global rdistnn_distance_engine "`smv_distance'"
global rdistnn_jackknife_engine "`smv_jackknife'"

// Canonical report arrays and route-specific retained controls belong to the
// current snapshot too. Corruption must fail before a consumer repairs it.
test_section "current retained-state reporting and route contracts"
clear
set obs 12
gen long smvc_id=_n
gen double smvc_x=(_n-6.5)/3
gen double smvc_z=cond(mod(_n,2),-.5,.5)
gen double smvc_y=1+smvc_x^2+sin(_n)/10
gen double smvc_den=.1+.01*smvc_x
gen byte smvc_d=mod(_n,2)
gen double smvc_e=cond(mod(_n,2),0,1)
gen double smvc_m0=.2+smvc_x
gen double smvc_m1=.5+smvc_x
matrix SMVC_TARGETS=(-.5,0\0,0\-.5,0)
foreach family in np cate clate {
    local extra ""
    local cases "target_identity target_dimensions estimates_dimensions ses_dimensions negative_se vcov_dimensions diagnostics_dimensions missing_estimate missing_se inconsistent_se negative_variance asymmetric_vcov missing_variance se_method jackknife_method"
    if "`family'"=="cate" {
        local extra "tcol(smvc_d) ecol(smvc_e) m0col(smvc_m0) m1col(smvc_m1)"
        local cases "bad_arm_link bad_propensity_link missing_arm_link missing_propensity_link negative_clip missing_clip large_clip"
    }
    if "`family'"=="clate" {
        local extra "pseudooutcome denominator(smvc_den)"
        local cases "negative_tolerance missing_tolerance status_dimensions components_dimensions componentcov_dimensions status_values status_missing status_availability inconsistent_se negative_variance asymmetric_vcov se_method jackknife_method status_ok_missing status_ok_zero status_zero_nonzero status_failed_finite"
    }
    quietly `family'_dnn smvc_y if !inlist(smvc_id,1,5,12), ///
        xcols(smvc_x smvc_z) x0matrix(SMVC_TARGETS) ///
        tuningmode(grid) sgrid(2 3 2) sdefault(3) `extra'
    estimates store SMVC_VALID
    foreach case of local cases {
        quietly estimates restore SMVC_VALID
        _smv_corrupt, case(`case')
        estimates store SMVC_INVALID
        preserve
        clear
        set obs 3
        gen str12 caller="preserved"
        gen double smvc_x=777
        quietly estimates restore SMVC_INVALID
        mata: mata clear
        _rdistnn_test_e_snapshot, name(SMVC_STATE)
        _smv_check_routes, estimates(SMVC_INVALID) snapshot(SMVC_STATE) case(`case') ///
            routes("replay results vcov diagnostics plot weights predict candidate") family(`family') query
        _rdistnn_test_e_snapshot_drop, name(SMVC_STATE)
        restore
        estimates drop SMVC_INVALID
    }
    estimates drop SMVC_VALID
}

test_section "complete current candidate grids validate before slicing"
foreach family in np cate clate {
    local extra ""
    if "`family'"=="cate" local extra "pseudooutcome"
    if "`family'"=="clate" local extra "pseudooutcome denominator(smvc_den)"
    quietly `family'_dnn smvc_y if !inlist(smvc_id,1,5,12), ///
        xcols(smvc_x smvc_z) x0matrix(SMVC_TARGETS) ///
        tuningmode(grid) sgrid(2 3 2 9) sdefault(3) `extra'
    estimates store SMVC_GRID_VALID
    local cases "grid_candidates_shape grid_estimates_shape grid_se_shape grid_diag_shape grid_default grid_id grid_computed grid_feasible grid_alias grid_availability grid_target_names grid_count"
    if "`family'"=="clate" local cases "`cases' grid_status_shape grid_y_shape grid_d_shape"
    foreach case of local cases {
        quietly estimates restore SMVC_GRID_VALID
        _smv_corrupt, case(`case')
        estimates store SMVC_GRID_INVALID
        _rdistnn_test_e_snapshot, name(SMVC_GRID_STATE)
        _smv_check_routes, estimates(SMVC_GRID_INVALID) snapshot(SMVC_GRID_STATE) case(`case') ///
            routes("replay results diagnostics plot predict candidate") family(`family') ///
            callertest("_N==12 & smvc_y[1]<.")
        _rdistnn_test_e_snapshot_drop, name(SMVC_GRID_STATE)
        estimates drop SMVC_GRID_INVALID
    }
    quietly estimates restore SMVC_GRID_VALID
    quietly rdistnn_candidate, candidate(3)
    local selected_candidate=r(candidate)
    local selected_order=r(s)
    test_assert "`family': valid aliased candidates remain reusable" `selected_candidate'==3 & `selected_order'==2
    estimates drop SMVC_GRID_VALID
    quietly `family'_dnn smvc_y, xcols(smvc_x smvc_z) x0matrix(SMVC_TARGETS) ///
        tuningmode(grid) lambdagrid(1e-20 .4 8e307) `extra'
    test_assert "`family': numerical infeasible optional orders retain current grid status" e(tuning_n_skipped)==2 & e(tuning_n_feasible)==1
}

test_section "paired component arrays validate at their shared retained-state owner"
foreach case in components_missing components_zero_denominator component_sign_flag ///
    componentcov_negative componentcov_asymmetric componentcov_partial ///
    componentcov_cross_missing componentcov_false_success componentcov_nose_finite {
    local seoption ""
    if "`case'"=="componentcov_nose_finite" local seoption "nose"
    quietly clate_dnn smvc_y, pseudooutcome denominator(smvc_den) ///
        xcols(smvc_x smvc_z) x0matrix(SMVC_TARGETS) s(3) `seoption'
    _smv_corrupt, case(`case')
    estimates store SMVC_COMPONENT_INVALID
    _rdistnn_test_e_snapshot, name(SMVC_COMPONENT_STATE)
    _smv_check_routes, estimates(SMVC_COMPONENT_INVALID) snapshot(SMVC_COMPONENT_STATE) ///
        case(`case') routes("components predict") family(clate) ///
        callertest("_N==12 & smvc_y[1]<.")
    _rdistnn_test_e_snapshot_drop, name(SMVC_COMPONENT_STATE)
    estimates drop SMVC_COMPONENT_INVALID
}
program drop _smv_corrupt
program drop _smv_check_routes

test_section "valid retained ratio uncertainty states"
clear
set obs 8
gen double state_x=_n
gen double state_den=1+_n/100
gen double state_y=2*state_den
quietly clate_dnn state_y, pseudooutcome denominator(state_den) xcols(state_x) x0(0) s(3)
capture quietly rdistnn_results
local rc=_rc
test_assert "defined zero-score ratio retains exact zero covariance status" `rc'==0 & e(status)[1,4]==2 & e(estimate_vcov)[1,1]==0
quietly clate_dnn state_y, pseudooutcome denominator(state_den) xcols(state_x) x0(0) s(3) nose
capture quietly rdistnn_results
local rc=_rc
test_assert "nose ratio retains intentionally uncomputed uncertainty" `rc'==0 & e(status)[1,4]==1 & missing(e(se_targets)[1,1])
replace state_den=0
replace state_y=_n^2
quietly clate_dnn state_y, pseudooutcome denominator(state_den) xcols(state_x) x0(0) s(3)
capture quietly rdistnn_results
local rc=_rc
test_assert "undefined ratio retains component-success status and missing ratio covariance" `rc'==0 & e(status)[1,1]==1 & e(status)[1,4]==0 & missing(e(estimate_vcov)[1,1])
replace state_den=1
replace state_y=cond(mod(_n,2),-1,1)*2^-540
quietly clate_dnn state_y, pseudooutcome denominator(state_den) xcols(state_x) x0(0) s(3)
capture quietly rdistnn_results
local rc=_rc
test_assert "numerical covariance failure retains a finite defined ratio" `rc'==0 & e(status)[1,1]==0 & e(status)[1,4]==3 & missing(e(estimate_vcov)[1,1])

replace state_y=2^900
replace state_den=2^-900
quietly clate_dnn state_y, pseudooutcome denominator(state_den) xcols(state_x) x0(0) s(3)
capture quietly rdistnn_results, type(components)
local rc=_rc
test_assert "finite constant components retain ratio-overflow and successful component covariance" `rc'==0 & e(status)[1,1]==2 & e(status)[1,4]==0 & !missing(e(components)[1,1],e(components)[1,2]) & e(component_vcov)[1,1]==0
replace state_y=(_n^2)*1e200
replace state_den=1
quietly clate_dnn state_y, pseudooutcome denominator(state_den) xcols(state_x) x0(0) s(3)
capture quietly rdistnn_results, type(components)
local rc=_rc
test_assert "component covariance range failure retains finite component means" `rc'==0 & e(status)[1,1]==0 & e(status)[1,4]==3 & !missing(e(components)[1,1],e(components)[1,2]) & missing(e(component_vcov)[1,1],e(component_vcov)[2,2])
replace state_y=cond(mod(_n,2),-1,1)*2^-450
replace state_den=2^-1000
quietly clate_dnn state_y, pseudooutcome denominator(state_den) xcols(state_x) x0(0) s(3)
capture quietly rdistnn_results, type(components)
local rc=_rc
test_assert "ratio variance range failure preserves finite component covariance" `rc'==0 & e(status)[1,1]==0 & e(status)[1,4]==3 & !missing(e(component_vcov)[1,1],e(component_vcov)[2,2]) & missing(e(estimate_vcov)[1,1])
