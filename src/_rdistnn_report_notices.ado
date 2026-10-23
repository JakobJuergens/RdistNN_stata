version 17
// One public reporting boundary reads already validated fit facts. Numerical
// kernels and quiet nested consumers do not own notice frequency.
program define _rdistnn_report_notices
    version 17
    syntax [, TARGETS(numlist integer >=1)]
    if !c(noisily) exit
    local m=e(m)
    local targets : list uniq targets
    if "`targets'"=="" mata: rdistnn_distance_warning(st_matrix("e(model_groups)"))
    else mata: rdistnn_distance_warning(st_matrix("e(model_groups)")[,strtoreal(tokens("`targets'"))])
    local concentration_column=cond("`e(cmd)'"=="np_tdnn",5,1)
    // The established near-uniform guide is descriptive, not a validity test.
    local count 0
    local examples ""
    forvalues j=1/`m' {
        if "`targets'"!="" {
            local included : list j in targets
            if !`included' continue
        }
        if e(N)>=2 & el(e(diagnostics),`j',`concentration_column')>=.9*e(N) {
            local ++count
            if `count'<=5 local examples "`examples' `j'"
        }
    }
    if `count' {
        local omitted ""
        if `count'>5 local omitted " (`=`count'-5' more)"
        local problem "weak localization: ESS is at least .9*n, the near-uniform-weight guide"
        local hint "Inspect e(diagnostics) and the target profiles; consider smaller k() or larger s() within feasible bounds. This guide is not a validity test."
        if "`e(cmd)'"=="np_tdnn" {
            local problem "weak localization: signed-weight concentration is at least .9*n"
            local hint "Inspect e(diagnostics) and the two localization orders. Inverse squared-weight concentration for signed TDNN weights is not probability-weight ESS or a validity test."
        }
        _rdistnn_messages warning, message("`count' target(s): `problem'. Examples:`examples'`omitted'.") hint("`hint'")
    }
    if "`e(tuning_mode)'"=="grid" {
        local count 0
        local examples ""
        forvalues g=1/`=rowsof(e(tuning_candidates))' {
            if el(e(tuning_candidates),`g',3)==0 & el(e(tuning_candidates),`g',4)==0 {
                local ++count
                if `count'<=5 local examples "`examples' `g'"
            }
        }
        if `count' {
            local omitted ""
            if `count'>5 local omitted " (`=`count'-5' more)"
            _rdistnn_messages warning, message("`count' infeasible optional tuning candidate(s) retained. Examples:`examples'`omitted'.") ///
                hint("Inspect e(tuning_candidates) and e(tuning_reason_codes); choose feasible orders with k <= s < n when requesting inference.")
        }
    }
    if "`e(cmd)'"=="clate_dnn" {
        forvalues family=1/4 {
            local count 0
            local examples ""
            forvalues j=1/`m' {
                if "`targets'"!="" {
                    local included : list j in targets
                    if !`included' continue
                }
                local affected=0
                if `family'==1 local affected=(el(e(status),`j',1)!=0)
                if `family'==2 local affected=(el(e(status),`j',2)==1)
                if `family'==3 local affected=(el(e(status),`j',3)==1)
                if `family'==4 & "`e(se_method)'"!="none" local affected=(el(e(status),`j',5)==0)
                if `affected' {
                    local ++count
                    if `count'<=5 local examples "`examples' `j'"
                }
            }
            if `count' {
                local omitted ""
                if `count'>5 local omitted " (`=`count'-5' more)"
                local kind "warning"
                if `family'==1 local problem "ratios unavailable"
                if `family'==2 {
                    local problem "negative first stages; defined signed ratios are retained"
                    local kind "note"
                }
                if `family'==3 local problem "first stages near zero; ratios are not regularized"
                if `family'==4 local problem "Wald intervals unavailable"
                _rdistnn_messages `kind', message("`count' target(s): `problem'. Examples:`examples'`omitted'.") ///
                    hint("See e(status) and e(components) before interpreting these targets.")
            }
        }
    }
end
