*! version 0.1.0  2026-10-18
// Private scalar/grid plan. No score, covariance, estimation post or data access.
program define _rdistnn_tuning_plan, rclass
    version 17
    syntax, N(integer) KX(integer) K(integer) S(integer) LAMBDA(real) ///
        MODE(string) SDEFAULT(integer) [SGRID(numlist) LAMBDAGRID(numlist) noSE]
    local tuning_mode "`mode'"
    if `lambda' <= 0 | missing(`lambda') {
        di as error "lambda() must be a positive finite number."
        exit 198
    }
    if "`lambdagrid'" != "" & "`sgrid'" != "" {
        di as error "Specify either lambdagrid() or sgrid(), not both."
        exit 198
    }
    if "`tuning_mode'" == "default" & ///
       ("`lambdagrid'" != "" | "`sgrid'" != "" | `sdefault' != 0) {
        di as error "lambdagrid(), sgrid(), and sdefault() require tuningmode(grid)."
        exit 198
    }
    if "`tuning_mode'" == "grid" & `s' != 0 {
        di as error "s() cannot be combined with tuningmode(grid)."
        exit 198
    }
    if `sdefault' < 0 {
        di as error "sdefault() must be a positive integer when specified."
        exit 198
    }
    if `sdefault' > 0 & "`sgrid'" == "" {
        di as error "sdefault() requires sgrid()."
        exit 198
    }

    local paper_base = `n'^((`kx' + 1) / (`kx' + 2))
    local tuning_source "paper"
    local candidate_s ""
    local candidate_lambda ""
    local default_index 1

    if "`tuning_mode'" == "default" {
        if `s' == 0 local s = floor(`lambda' * `paper_base')
        else local tuning_source "s"
        local default_s = `s'
        local default_lambda = `lambda'
        local candidate_s "`s'"
        local candidate_lambda "`lambda'"
    }
    else if "`sgrid'" != "" {
        local tuning_source "sgrid"
        foreach sg of local sgrid {
            if `sg' < 1 | `sg' != floor(`sg') | missing(`sg') {
                di as error "sgrid() must contain positive integers."
                exit 198
            }
            local candidate_s "`candidate_s' `sg'"
            local candidate_lambda "`candidate_lambda' ."
        }
        local candidate_s : list retokenize candidate_s
        local candidate_lambda : list retokenize candidate_lambda
        local first_s : word 1 of `candidate_s'
        if `sdefault' == 0 local default_s = `first_s'
        else local default_s = `sdefault'
        local default_lambda = .
        local default_index 0
        local g = 0
        foreach sg of local candidate_s {
            local ++g
            if `default_index' == 0 & `sg' == `default_s' local default_index = `g'
        }
        if `default_index' == 0 {
            local candidate_s "`default_s' `candidate_s'"
            local candidate_lambda ". `candidate_lambda'"
            local default_index 1
        }
    }
    else {
        local tuning_source "lambda_grid"
        if "`lambdagrid'" == "" local lambdagrid ".3 .4 .5 .75"
        foreach lg of local lambdagrid {
            if `lg' <= 0 | missing(`lg') {
                di as error "lambdagrid() must contain positive finite numbers."
                exit 198
            }
            local sg = floor(`lg' * `paper_base')
            local candidate_s "`candidate_s' `sg'"
            local candidate_lambda "`candidate_lambda' `lg'"
        }
        local candidate_s : list retokenize candidate_s
        local candidate_lambda : list retokenize candidate_lambda
        local default_s = floor(`lambda' * `paper_base')
        local default_lambda = `lambda'
        local default_index 0
        local g = 0
        foreach sg of local candidate_s {
            local ++g
            local lg : word `g' of `candidate_lambda'
            if `default_index' == 0 & `sg' == `default_s' & ///
               abs(`lg' - `default_lambda') <= sqrt(c(epsdouble)) * max(1, abs(`default_lambda')) {
                local default_index = `g'
            }
        }
        if `default_index' == 0 {
            local g = 0
            foreach sg of local candidate_s {
                local ++g
                if `default_index' == 0 & `sg' == `default_s' local default_index = `g'
            }
        }
        if `default_index' == 0 {
            local candidate_s "`default_s' `candidate_s'"
            local candidate_lambda ". `candidate_lambda'"
            local default_index 1
        }
    }

    local s = `default_s'
    local ncandidates : word count `candidate_s'
    tempname candidate_mat
    matrix `candidate_mat' = J(`ncandidates', 8, .)
    local candidate_names ""
    local seen_s ""
    local nfeasible 0
    local nunique 0
    local nunique_feasible 0
    forvalues g = 1/`ncandidates' {
        local sg : word `g' of `candidate_s'
        local lg : word `g' of `candidate_lambda'
        local feasible 1
        local reason 0
        if `k' > `sg' {
            local feasible 0
            local reason 1
        }
        else if "`se'" == "" & `sg' >= `n' {
            local feasible 0
            local reason 2
        }
        else if "`se'" != "" & `sg' > `n' {
            local feasible 0
            local reason 3
        }
        local computed_id 0
        local h = 0
        foreach old_s of local seen_s {
            local ++h
            if `computed_id' == 0 & `old_s' == `sg' local computed_id = `h'
        }
        if `computed_id' == 0 {
            local ++nunique
            local computed_id = `nunique'
            local seen_s "`seen_s' `sg'"
            if `feasible' local ++nunique_feasible
        }
        local nfeasible = `nfeasible' + `feasible'
        matrix `candidate_mat'[`g', 1] = `g'
        matrix `candidate_mat'[`g', 2] = `computed_id'
        matrix `candidate_mat'[`g', 3] = (`g' == `default_index')
        matrix `candidate_mat'[`g', 4] = `feasible'
        matrix `candidate_mat'[`g', 5] = `reason'
        matrix `candidate_mat'[`g', 6] = `feasible'
        matrix `candidate_mat'[`g', 7] = `lg'
        matrix `candidate_mat'[`g', 8] = `sg'
        local candidate_names "`candidate_names' candidate`g'"
    }
    local candidate_names : list retokenize candidate_names
    matrix rownames `candidate_mat' = `candidate_names'
    matrix colnames `candidate_mat' = candidate_id computed_id is_default feasible reason status lambda s

    if el(`candidate_mat', `default_index', 4) == 0 {
        local default_reason = el(`candidate_mat', `default_index', 5)
        if `default_reason' == 1 {
            di as error "Default s(`s') < k(`k'). Increase n or lambda(), or decrease k()."
        }
        else if `default_reason' == 2 {
            di as error "Delete-one inference requires s() < the estimation sample size; got s = `s' and n = `n'. Reduce s() or lambda(), or specify nose."
        }
        else {
            di as error "Default s(`s') must be no greater than n(`n') when nose is specified."
        }
        exit 198
    }

    return scalar s=`s'
    return scalar default_lambda=`default_lambda'
    return scalar default_index=`default_index'
    return scalar ncandidates=`ncandidates'
    return scalar nfeasible=`nfeasible'
    return scalar nunique_feasible=`nunique_feasible'
    return local tuning_source "`tuning_source'"
    return local candidate_names "`candidate_names'"
    return matrix candidates=`candidate_mat'
end
