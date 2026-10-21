// Common diagnostics contract. These tests catch dropped/reordered target IDs,
// missing common returns, lost scientific fields and saved-data side effects.
version 17
do "tests/test_framework.do"

test_section "Estimator diagnostics: common selection, localization and saving"
clear
set obs 32
gen double x1 = (_n-16)/8
gen double x2 = sin(_n/3)
gen double signal = 1 + x1 + x2^2
gen double denominator = .8 + .1*x2
matrix DIAG_X0 = (0,.3 \ .5,-.5 \ -1,1)
local common_names "target_id x0_1 x0_2 ess max_weight top1_weight_share top5_weight_share top10_weight_share n_nonzero_weight distance_available radius_rank_k radius_rank_s weighted_mean_distance weighted_median_distance"

foreach family in cate clate np {
    if "`family'"=="cate" quietly cate_dnn signal, pseudooutcome xcols(x1 x2) x0matrix(DIAG_X0) tuningmode(grid)
    if "`family'"=="clate" quietly clate_dnn signal, denominator(denominator) pseudooutcome xcols(x1 x2) x0matrix(DIAG_X0) tuningmode(grid)
    if "`family'"=="np" quietly np_dnn signal, xcols(x1 x2) x0matrix(DIAG_X0) tuningmode(grid)
    matrix DIAG_RAW = e(diagnostics)
    matrix DIAG_TARGETS = e(targets)
    matrix DIAG_B = e(b)
    matrix DIAG_V = e(V)
    matrix DIAG_META = e(tuning_candidates)
    matrix DIAG_CANDRAW = e(tuning_diagnostics)
    local original_n = e(N)
    // estimates store retains an e(sample) marker in the dataset; establish
    // the data baseline after that deliberate setup change.
    estimates store diag_`family'
    quietly datasignature
    local original_data "`r(datasignature)'"

    capture noisily `family'_dnn_diagnostics
    local rc = _rc
    matrix DIAG_LOC = r(localization)
    matrix DIAG_FULL = r(diagnostics)
    matrix DIAG_COMPACT = r(compact)
    local scope "`r(scope)'"
    local nrows = r(n_rows)
    local locnames : colnames DIAG_LOC
    test_assert "[`family' diagnostics] default target localization has full common schema" `rc'==0 & rowsof(DIAG_LOC)==3 & colsof(DIAG_LOC)==14 & "`locnames'"=="`common_names'"
    test_assert "[`family' diagnostics] default table retains all targets" `nrows'==3 & "`scope'"=="targets" & rowsof(DIAG_FULL)==3

    tempfile diag_saved
    capture noisily `family'_dnn_diagnostics, targets(3 1) detail saving("`diag_saved'") replace
    local rc = _rc
    matrix DIAG_LOC = r(localization)
    matrix DIAG_FULL = r(diagnostics)
    matrix DIAG_COMPACT = r(compact)
    local scope "`r(scope)'"
    local nrows = r(n_rows)
    local selected "`r(selected)'"
    local locnames : colnames DIAG_LOC
    local locrows : rownames DIAG_LOC
    local fullrows : rownames DIAG_FULL
    if "`family'"=="clate" {
        matrix DIAG_COMPONENTS = r(components)
        matrix DIAG_STATUS = r(status)
    }
    test_assert "[`family' diagnostics] accepts reordered target selection" `rc'==0
    test_assert "[`family' diagnostics] selected localization has actual 2 by 14 dimensions" rowsof(DIAG_LOC)==2 & colsof(DIAG_LOC)==14 & "`locnames'"=="`common_names'"
    test_assert "[`family' diagnostics] IDs, coordinates and original order survive selection" DIAG_LOC[1,1]==3 & DIAG_LOC[2,1]==1 & DIAG_LOC[1,2]==-1 & DIAG_LOC[1,3]==1 & DIAG_LOC[2,2]==0 & DIAG_LOC[2,3]==.3
    test_assert "[`family' diagnostics] selected target row labels retain identities" "`locrows'"=="target3 target1" & "`fullrows'"=="target3 target1"
    test_assert "[`family' diagnostics] selected metadata reports original identities" `nrows'==2 & "`scope'"=="targets" & "`selected'"=="3 1"
    if `rc'==0 & rowsof(DIAG_LOC)==2 & colsof(DIAG_LOC)==14 {
        mata: st_numscalar("diag_common_error", max(abs(st_matrix("DIAG_LOC")[,4..14] - st_matrix("DIAG_RAW")[(3,1),1..11])))
        test_assert "[`family' diagnostics] localization preserves all eleven stored common fields" scalar(diag_common_error)<1e-12
        preserve
        use "`diag_saved'", clear
        local saved_n = _N
        local saved_ids = target_id[1]==3 & target_id[2]==1
        local saved_coords = x0_1[1]==-1 & x0_2[1]==1 & x0_1[2]==0 & x0_2[2]==.3
        local saved_extras 0
        if "`family'"=="cate" {
            capture confirm variable treated_mass control_mass treated_ess control_ess physical_radius_90 top5_score_share
            local saved_extras = _rc==0
        }
        if "`family'"=="clate" {
            capture confirm variable instrument_1_mass instrument_0_mass instrument_1_ess instrument_0_ess first_stage denominator_stat ratio_status covariance_status wald_available
            local saved_extras = _rc==0
        }
        if "`family'"=="np" {
            capture confirm variable weighted_mean_y weighted_sd_y weighted_y_q10 weighted_y_q50 weighted_y_q90
            local saved_extras = _rc==0
        }
        restore
        test_assert "[`family' diagnostics] saving retains selected identity, coordinates and family fields" `saved_n'==2 & `saved_ids' & `saved_coords' & `saved_extras'
    }
    else test_assert "[`family' diagnostics] selected saving returns usable common localization" 0
    local fullcols = cond("`family'"=="np",20,28)
    local compactcols = cond("`family'"=="np",6,7)
    test_assert "[`family' diagnostics] legacy full/compact tables retain fields and IDs" rowsof(DIAG_FULL)==2 & colsof(DIAG_FULL)==`fullcols' & DIAG_FULL[1,1]==3 & DIAG_FULL[2,1]==1 & rowsof(DIAG_COMPACT)==2 & colsof(DIAG_COMPACT)==`compactcols'
    if "`family'"=="clate" {
        local fullnames : colnames DIAG_FULL
        test_assert "[clate diagnostics] selected components and statuses retain full schema" rowsof(DIAG_COMPONENTS)==2 & colsof(DIAG_COMPONENTS)==9 & DIAG_COMPONENTS[1,1]==3 & DIAG_COMPONENTS[2,1]==1 & rowsof(DIAG_STATUS)==2 & colsof(DIAG_STATUS)==5 & strpos("`fullnames'","instrument_1_mass")>0
        matrix DIAG_ORIG_STATUS = e(status)
        matrix DIAG_ORIG_COMPONENTS = e(components)
        test_assert "[clate diagnostics] scientific tables use original target indices" DIAG_COMPONENTS[1,2]==DIAG_ORIG_COMPONENTS[3,1] & DIAG_COMPONENTS[2,3]==DIAG_ORIG_COMPONENTS[1,2] & DIAG_STATUS[1,1]==DIAG_ORIG_STATUS[3,1]
    }

    capture noisily `family'_dnn_diagnostics, candidates(3 1)
    local rc = _rc
    matrix DIAG_CAND = r(diagnostics)
    matrix DIAG_CANDCOMPACT = r(compact)
    local scope "`r(scope)'"
    local selected "`r(selected)'"
    local nrows = r(n_rows)
    test_assert "[`family' diagnostics] candidates retain legacy schemas and requested order" `rc'==0 & rowsof(DIAG_CAND)==2 & colsof(DIAG_CAND)==13 & rowsof(DIAG_CANDCOMPACT)==2 & colsof(DIAG_CANDCOMPACT)==10 & DIAG_CAND[1,1]==3 & DIAG_CAND[2,1]==1 & "`scope'"=="candidates" & `nrows'==2 & "`selected'"=="3 1"
    if `rc'==0 {
        mata: st_numscalar("diag_candidate_error", max(abs(st_matrix("DIAG_CAND") - (st_matrix("DIAG_META")[(3,1),], st_matrix("DIAG_CANDRAW")[(3,1),]))))
        test_assert "[`family' diagnostics] selected candidate metadata and diagnostics unchanged" scalar(diag_candidate_error)<1e-12
    }
    tempfile diag_candidates
    capture noisily `family'_dnn_diagnostics, allcandidates saving("`diag_candidates'") replace
    local rc = _rc
    matrix DIAG_ALL = r(diagnostics)
    local nrows = r(n_rows)
    test_assert "[`family' diagnostics] allcandidates keeps every stored row" `rc'==0 & `nrows'==rowsof(DIAG_META) & rowsof(DIAG_ALL)==rowsof(DIAG_META) & colsof(DIAG_ALL)==13
    if `rc'==0 {
        preserve
        use "`diag_candidates'", clear
        local candidate_save_ok = _N==rowsof(DIAG_META) & candidate_id[1]==1 & candidate_id[_N]==rowsof(DIAG_META)
        capture confirm variable computed_id feasible reason status lambda s ess max_weight top_n top5_weight_share degenerate
        local candidate_save_ok = `candidate_save_ok' & _rc==0
        restore
        test_assert "[`family' diagnostics] candidate saving preserves full fields and identities" `candidate_save_ok'
    }

    capture noisily `family'_dnn_diagnostics, targets(1) candidates(1)
    test_assert "[`family' diagnostics] target/candidate selections are mutually exclusive" _rc==198
    capture noisily `family'_dnn_diagnostics, targets(1) allcandidates
    test_assert "[`family' diagnostics] target/allcandidate selections are mutually exclusive" _rc==198
    capture noisily `family'_dnn_diagnostics, candidates(1) allcandidates
    test_assert "[`family' diagnostics] candidate selectors are mutually exclusive" _rc==198
    capture noisily `family'_dnn_diagnostics, targets(4)
    test_assert "[`family' diagnostics] target range checked" _rc==198
    capture noisily `family'_dnn_diagnostics, candidates(999)
    test_assert "[`family' diagnostics] candidate range checked" _rc==198
    capture noisily `family'_dnn_diagnostics, saving("`diag_candidates'")
    test_assert "[`family' diagnostics] saving refuses overwrite without replace" _rc==602
    quietly datasignature
    local final_data "`r(datasignature)'"
    matrix DIAG_AFTER_RAW = e(diagnostics)
    matrix DIAG_AFTER_B = e(b)
    matrix DIAG_AFTER_V = e(V)
    mata: st_numscalar("diag_estimation_unchanged", (st_matrix("DIAG_AFTER_RAW")==st_matrix("DIAG_RAW")) & (st_matrix("DIAG_AFTER_B")==st_matrix("DIAG_B")) & (st_matrix("DIAG_AFTER_V")==st_matrix("DIAG_V")))
    test_assert "[`family' diagnostics] successful and rejected calls preserve e() and data" scalar(diag_estimation_unchanged) & "`e(cmd)'"=="`family'_dnn" & e(N)==`original_n' & "`original_data'"=="`final_data'"
}

// Loading a restored result must not depend on which public command ran last.
estimates restore diag_cate
capture noisily cate_dnn_diagnostics, targets(2)
local restore_rc = _rc
matrix DIAG_RESTORED = r(localization)
test_assert "[diagnostics] stored results restore target localization" `restore_rc'==0 & rowsof(DIAG_RESTORED)==1 & colsof(DIAG_RESTORED)==14 & DIAG_RESTORED[1,1]==2

// Undefined CLATE ratios still retain their localization and status rows.
gen double zero_denominator = 0
quietly clate_dnn signal, denominator(zero_denominator) pseudooutcome xcols(x1 x2) x0matrix(DIAG_X0) s(8) nose
capture noisily clate_dnn_diagnostics, targets(3 1)
local undefined_rc = _rc
matrix DIAG_UNDEFINED = r(localization)
matrix DIAG_UNDEFINED_STATUS = r(status)
matrix DIAG_UNDEFINED_COMPONENTS = r(components)
test_assert "[clate diagnostics] undefined ratios keep canonical target/status identity" `undefined_rc'==0 & rowsof(DIAG_UNDEFINED)==2 & colsof(DIAG_UNDEFINED)==14 & DIAG_UNDEFINED[1,1]==3 & rowsof(DIAG_UNDEFINED_STATUS)==2 & colsof(DIAG_UNDEFINED_STATUS)==5 & DIAG_UNDEFINED_STATUS[1,1]==1 & DIAG_UNDEFINED_COMPONENTS[1,1]==3 & missing(DIAG_UNDEFINED_COMPONENTS[1,9])
