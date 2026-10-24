// Batch, tuning, diagnostics, plotting, and deterministic R parity for np_dnn.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
quietly do "tests/fixtures/r_priority_fixtures.do"

clear
input double x1 double x2
-1     0
 0    -1
 0     0
 1     0
 0     1
 1     1
-1    -1
 2    -0.5
-0.5   1.5
 0.75 -1.25
 1.5   0.5
-1.5   0.75
end
gen long fixture_i = _n
gen double y_np = 0.7*x1 - 0.4*x2 + fixture_i/20
matrix X0_np = (0.2,-0.1 \ -0.5,0.75 \ 1.2,0.4)

test_section "np_dnn: [R parity] batch targets and diagnostics"

quietly np_dnn y_np, xcols(x1 x2) x0matrix(X0_np) k(2) s(5)
matrix B_np_batch = e(b)
matrix V_np_batch = e(V)
matrix E_np_batch = e(ess_targets)
matrix D_np_batch = e(diagnostics)
matrix T_np_batch = e(targets)
mata {
    st_numscalar("_np_batch_mu_r", max(abs(st_matrix("B_np_batch") - st_matrix("R_np_batch_exact_mu")')))
    st_numscalar("_np_batch_se_r", max(abs(sqrt(diagonal(st_matrix("V_np_batch"))) - st_matrix("R_np_batch_exact_se"))))
    st_numscalar("_np_batch_ess_r", max(abs(st_matrix("E_np_batch") - st_matrix("R_np_batch_exact_ess")')))
    st_numscalar("_np_batch_diag_r", max(abs(st_matrix("D_np_batch") - st_matrix("R_np_batch_exact_diagnostics"))))
    st_numscalar("_np_batch_median_distance_r", max(abs(
        st_matrix("D_np_batch")[,11] - st_matrix("R_np_batch_exact_diagnostics")[,11])))
    st_numscalar("_np_batch_outcome_quantiles_r", max(abs(
        st_matrix("D_np_batch")[,15..17] - st_matrix("R_np_batch_exact_diagnostics")[,15..17])))
    st_numscalar("_np_batch_vsym", max(abs(st_matrix("V_np_batch") - st_matrix("V_np_batch")')))
}
test_assert "[R parity] batch estimates match RdistNN" scalar(_np_batch_mu_r) < 1e-11
test_assert "[R parity] batch standard errors match RdistNN" scalar(_np_batch_se_r) < 1e-11
test_assert "[R parity] batch ESS matches RdistNN" scalar(_np_batch_ess_r) < 1e-11
test_assert "[R parity] target diagnostics match RdistNN" scalar(_np_batch_diag_r) < 1e-11
test_assert "[R parity] weighted median distances match RdistNN" ///
    scalar(_np_batch_median_distance_r) < 1e-11
test_assert "[R parity] NP outcome quantiles match RdistNN" ///
    scalar(_np_batch_outcome_quantiles_r) < 1e-11
test_assert "[batch] e(V) is symmetric" scalar(_np_batch_vsym) < 1e-12
test_assert "[batch] e(b) has one coefficient per target" colsof(B_np_batch) == 3
test_assert "[batch] e(V) retains cross-target covariance" rowsof(V_np_batch) == 3 & colsof(V_np_batch) == 3
test_assert "[batch] targets preserve input row order" T_np_batch[2,1] == -0.5 & T_np_batch[3,2] == 0.4
test_assert "[naming] xcols() is recorded as x_cols" "`e(x_cols)'" == "x1 x2"
test_assert "[naming] x0matrix() records R X0 source" "`e(target_source)'" == "X0"
test_assert "[meta] batch target count is recorded" e(m) == 3

matrix X0_np_reordered = (1.2,0.4 \ 0.2,-0.1 \ -0.5,0.75)
quietly np_dnn y_np, xcols(x1 x2) x0matrix(X0_np_reordered) k(2) s(5)
matrix B_np_reordered = e(b)
matrix V_np_reordered = e(V)
mata {
    _np_perm = (3,1,2)
    st_numscalar("_np_batch_perm_b", max(abs(st_matrix("B_np_reordered") - st_matrix("B_np_batch")[1,_np_perm])))
    st_numscalar("_np_batch_perm_v", max(abs(st_matrix("V_np_reordered") - st_matrix("V_np_batch")[_np_perm,_np_perm])))
}
test_assert "[determinism] target reordering permutes estimates only" scalar(_np_batch_perm_b) < 1e-12
test_assert "[determinism] target reordering permutes joint covariance" scalar(_np_batch_perm_v) < 1e-12


test_section "np_dnn: [nose] estimate-only batch behavior"

quietly np_dnn y_np, xcols(x1 x2) x0matrix(X0_np) k(2) s(12) nose
matrix B_np_nose = e(b)
local nose_method "`e(se_method)'"
capture matrix list e(V)
local nose_v_rc = _rc
test_assert "[nose] s=n is feasible without jackknife" e(s) == e(N)
test_assert "[nose] no covariance matrix is posted" `nose_v_rc' != 0
test_assert "[nose] no-SE method is recorded" "`nose_method'" == "none"
test_assert "[nose] batch estimates remain finite" B_np_nose[1,1] < . & B_np_nose[1,3] < .
local nose_b_before = B_np_nose[1,1]
quietly np_dnn
matrix B_np_nose_replay = e(b)
test_assert "[nose] replay preserves estimate-only results" ///
    abs(B_np_nose_replay[1,1] - `nose_b_before') < 1e-12 & ///
    "`e(se_method)'" == "none" & "`e(jackknife)'" == "none"
quietly np_dnn_plot, type(estimates) name(np_test_nose) replace nodraw
test_assert "[nose] estimate plot works without e(V)" r(n_rows) == 3

quietly np_dnn y_np, xcols(x1 x2) x0(0.2 -0.1) k(2) nose ///
    tuningmode(grid) sgrid(6 12) sdefault(12)
matrix C_np_nose_grid = e(tuning_candidates)
capture matrix list e(V)
local nose_grid_v_rc = _rc
test_assert "[nose] grid permits s=n candidate" C_np_nose_grid[2,4] == 1 & e(s) == e(N)
test_assert "[nose] grid at s=n still omits e(V)" `nose_grid_v_rc' != 0

quietly np_dnn y_np, xcols(x1 x2) x0(0.2 -0.1) k(2) s(5)
matrix B_np_with_se = e(b)
quietly np_dnn y_np, xcols(x1 x2) x0(0.2 -0.1) k(2) s(5) nose
matrix B_np_without_se = e(b)
test_assert "[nose] estimate equals ordinary fit at the same k and s" ///
    abs(B_np_with_se[1,1] - B_np_without_se[1,1]) < 1e-12

test_section "np_dnn: [R parity] integrated paper tuning grid"

quietly np_dnn y_np, xcols(x1 x2) x0matrix(X0_np) k(2) tuningmode(grid)
matrix C_np_grid = e(tuning_candidates)
matrix G_np_grid = e(tuning_estimates)
matrix S_np_grid = e(tuning_se)
matrix D_np_grid = e(tuning_diagnostics)
matrix B_np_grid = e(b)
matrix V_np_grid = e(V)
matrix ESS_np_grid = e(ess_targets)
mata {
    _NEg = st_matrix("G_np_grid"); _NEr = st_matrix("R_np_grid_estimates")
    _NSg = st_matrix("S_np_grid"); _NSr = st_matrix("R_np_grid_se")
    _NDg = st_matrix("D_np_grid"); _NDr = st_matrix("R_np_grid_diagnostics")
    st_numscalar("_np_grid_c_r", max(abs(st_matrix("C_np_grid") - st_matrix("R_np_grid_candidates"))))
    st_numscalar("_np_grid_e_r", max(abs(editmissing(_NEg,0) - editmissing(_NEr,0))))
    st_numscalar("_np_grid_e_missing", sum(vec((_NEg:>=.) :!= (_NEr:>=.))))
    st_numscalar("_np_grid_s_r", max(abs(editmissing(_NSg,0) - editmissing(_NSr,0))))
    st_numscalar("_np_grid_s_missing", sum(vec((_NSg:>=.) :!= (_NSr:>=.))))
    st_numscalar("_np_grid_d_r", max(abs(editmissing(_NDg,0) - editmissing(_NDr,0))))
    st_numscalar("_np_grid_d_missing", sum(vec((_NDg:>=.) :!= (_NDr:>=.))))
    st_numscalar("_np_grid_default_b", max(abs(st_matrix("B_np_grid") - st_matrix("R_np_grid_estimates")[2,.])))
    st_numscalar("_np_grid_default_s", max(abs(sqrt(diagonal(st_matrix("V_np_grid")))' - st_matrix("R_np_grid_se")[2,.])))
    st_numscalar("_np_grid_default_ess", max(abs(st_matrix("ESS_np_grid") - st_matrix("R_np_grid_default_ess")')))
}
test_assert "[R parity] NP candidate order, lambda, and resolved s match RdistNN" scalar(_np_grid_c_r) < 1e-12
test_assert "[R parity] all NP candidate estimates match RdistNN" scalar(_np_grid_e_r) < 1e-11 & scalar(_np_grid_e_missing) == 0
test_assert "[R parity] all NP candidate SEs match RdistNN" scalar(_np_grid_s_r) < 1e-11 & scalar(_np_grid_s_missing) == 0
test_assert "[R parity] all NP candidate diagnostics match RdistNN" scalar(_np_grid_d_r) < 1e-11 & scalar(_np_grid_d_missing) == 0
test_assert "[R parity] NP default target ESS matches RdistNN" scalar(_np_grid_default_ess) < 1e-11
test_assert "[grid] default estimate is candidate 2" scalar(_np_grid_default_b) < 1e-11
test_assert "[grid] default variance diagonal is candidate 2" scalar(_np_grid_default_s) < 1e-11
test_assert "[grid] finite-sample default is lambda=0.4 and s=2" abs(e(lambda)-.4) < 1e-12 & e(s) == 2
test_assert "[grid] lambda=0.3 candidate is infeasible for k=2" C_np_grid[1,4] == 0 & C_np_grid[1,5] == 1
test_assert "[grid] three candidates are feasible and unique" e(tuning_n_feasible) == 3 & e(tuning_n_unique) == 3
test_assert "[grid] mode is recorded" "`e(tuning_mode)'" == "grid"

forvalues g = 2/4 {
    local sg = C_np_grid[`g',8]
    quietly np_dnn y_np, xcols(x1 x2) x0matrix(X0_np) k(2) s(`sg')
    matrix B_np_scalar = e(b)
    mata: st_numscalar("_np_grid_scalar", max(abs(st_matrix("B_np_scalar") - st_matrix("R_np_grid_estimates")[`g',.])))
    test_assert "[grid] candidate `g' equals scalar fit" scalar(_np_grid_scalar) < 1e-11
}

quietly np_dnn y_np, xcols(x1 x2) x0(0.2 -0.1) k(2) ///
    tuningmode(grid) sgrid(4 6 6 8) sdefault(6)
matrix C_np_dup = e(tuning_candidates)
matrix E_np_dup = e(tuning_estimates)
test_assert "[grid] duplicate s values share computed_id" C_np_dup[2,2] == C_np_dup[3,2]
test_assert "[grid] duplicate candidate estimates agree" abs(E_np_dup[2,1] - E_np_dup[3,1]) < 1e-12
test_assert "[grid] duplicate orders count once" e(tuning_n_unique) == 3

quietly np_dnn y_np, xcols(x1 x2) x0(0.2 -0.1) k(1) lambda(1)
test_assert "[override] NP lambda(1) retains the aggressive benchmark" abs(e(lambda)-1) < 1e-12 & e(s) == 6
quietly np_dnn y_np, xcols(x1 x2) x0(0.2 -0.1) k(1) ///
    tuningmode(grid) lambda(1) lambdagrid(.4 1)
matrix C_np_lambda_override = e(tuning_candidates)
test_assert "[override] NP custom lambdagrid() order is preserved" C_np_lambda_override[1,7] == .4 & C_np_lambda_override[2,7] == 1
test_assert "[override] NP lambda(1) selects its custom-grid row" e(tuning_default_candidate) == 2 & e(s) == 6

quietly np_dnn y_np, xcols(x1 x2) x0(0.2 -0.1) k(2) ///
    tuningmode(grid) sgrid(1 6 12) sdefault(6)
matrix C_np_skip = e(tuning_candidates)
test_assert "[grid] k>s candidate is retained as infeasible" C_np_skip[1,4] == 0 & C_np_skip[1,5] == 1
test_assert "[grid] jackknife s=n candidate is infeasible" C_np_skip[3,4] == 0 & C_np_skip[3,5] == 2
test_assert "[grid] skipped candidate count is recorded" e(tuning_n_skipped) == 2

test_section "np_dnn: [postestimation] diagnostics and plots"

quietly np_dnn y_np, xcols(x1 x2) x0matrix(X0_np) k(2) tuningmode(grid)
quietly np_dnn_diagnostics
matrix NP_diag_table = r(diagnostics)
matrix NP_diag_compact = r(compact)
local np_diag_scope "`r(scope)'"
test_assert "[diagnostics] target table has id, targets, and 17 R fields" colsof(NP_diag_table) == 20
test_assert "[diagnostics] compact target table has six fields" colsof(NP_diag_compact) == 6
test_assert "[diagnostics] target row count matches e(m)" rowsof(NP_diag_table) == 3
test_assert "[diagnostics] default scope is targets" "`np_diag_scope'" == "targets"
test_assert "[diagnostics] weighted mean equals fitted estimate" abs(NP_diag_table[1,16] - R_np_grid_estimates[2,1]) < 1e-11

quietly np_dnn_diagnostics, allcandidates
matrix NP_candidate_diag = r(diagnostics)
local np_candidate_scope "`r(scope)'"
test_assert "[diagnostics] candidate table uses 13-field contract" colsof(NP_candidate_diag) == 13
test_assert "[diagnostics] all four candidates are returned" rowsof(NP_candidate_diag) == 4
test_assert "[diagnostics] candidate scope is recorded" "`np_candidate_scope'" == "candidates"

tempfile np_diag_saved
quietly np_dnn_diagnostics, saving(`np_diag_saved') replace
preserve
quietly use `np_diag_saved', clear
quietly count
local np_saved_n = r(N)
capture confirm variable weighted_y_q90
local np_saved_var_rc = _rc
restore
test_assert "[diagnostics] saving() writes one row per target" `np_saved_n' == 3
test_assert "[diagnostics] saving() preserves outcome quantiles" `np_saved_var_rc' == 0

foreach p in estimates diagnostics tuning weights outcome dashboard {
    local extra ""
    if "`p'" == "weights" local extra "view(kernel)"
    quietly np_dnn_plot, type(`p') `extra' name(np_test_`p') replace nodraw
    local plotted = r(n_rows)
    if "`p'"=="dashboard" local plotted = r(n_panels)
    local plot_type "`r(type)'"
    test_assert "[plot] type(`p') builds without drawing" `plotted' > 0 & "`plot_type'" == "`p'"
}

quietly np_dnn_plot, type(estimates) allcandidates name(np_test_est_candidates) replace nodraw
test_assert "[plot] candidate estimate paths are supported" r(n_rows) == 9
quietly np_dnn_plot, type(diagnostics) allcandidates name(np_test_diag_candidates) replace nodraw
test_assert "[plot] candidate diagnostic paths are supported" r(n_rows) == 3
quietly np_dnn_plot, type(weights) allcandidates name(np_test_weight_candidates) replace nodraw
test_assert "[plot] candidate weight profiles are supported" r(n_rows) == 3*e(N)*e(m)

test_section "np_dnn: [errors] expanded interface validation"

capture noisily np_dnn y_np, xcols(x1 x2) x0(0 0) x0matrix(X0_np)
test_assert "[errors] x0() and x0matrix() are mutually exclusive" _rc == 198
matrix X0_bad_np = (0,0,0)
capture noisily np_dnn y_np, xcols(x1 x2) x0matrix(X0_bad_np)
test_assert "[errors] target column mismatch is rejected" _rc == 198
matrix X0_missing_np = (0,. )
capture noisily np_dnn y_np, xcols(x1 x2) x0matrix(X0_missing_np)
test_assert "[errors] missing target values are rejected" _rc == 198
capture noisily np_dnn y_np, xcols(x1 x2) x0(0 0) s(5) tuningmode(grid)
test_assert "[errors] grid mode rejects explicit s()" _rc == 198
capture noisily np_dnn y_np, xcols(x1 x2) x0(0 0) lambdagrid(.3 .4)
test_assert "[errors] lambdagrid() requires grid mode" _rc == 198
capture noisily np_dnn y_np, xcols(x1 x2) x0(0 0) ///
    tuningmode(grid) lambdagrid(.3 .4) sgrid(4 6)
test_assert "[errors] lambda and s grids are mutually exclusive" _rc == 198

quietly regress y_np x1 x2
capture noisily np_dnn_diagnostics
test_assert "[errors] diagnostics require active np_dnn estimates" _rc == 301
capture noisily np_dnn_plot
test_assert "[errors] plots require active np_dnn estimates" _rc == 301
