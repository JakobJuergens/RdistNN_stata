// test-cate-dnn-postestimation.do
// Diagnostics tables and Stata-native CATE plot builders.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "cate_dnn postestimation: diagnostics tables"

rdnn_dgp_basic, n(120) seed(2310)
matrix X0_post = (-1,0 \ -.5,0 \ 0,0 \ .5,0 \ 1,0)
quietly cate_dnn y, xcols(x1 x2) x0matrix(X0_post) tcol(D) ecol(e) ///
    m0col(m0) m1col(m1) k(2) tuningmode(grid)
matrix B_post = e(b)
matrix V_post = e(V)
local s_post = e(s)
mata {
    _post_plot_data = rdistnn_plot_est_default(
        st_matrix("X0_post"), st_matrix("B_post")',
        sqrt(diagonal(st_matrix("V_post"))), `s_post',
        invnormal(.975))
    st_numscalar("_post_ci_width", min(_post_plot_data[,cols(_post_plot_data)] :-
        _post_plot_data[,cols(_post_plot_data)-1]))
}
test_assert "[plot data] stored SEs produce positive confidence-band width" scalar(_post_ci_width) > 0

quietly cate_dnn_diagnostics
matrix D_post = r(diagnostics)
local diag_scope "`r(scope)'"
local diag_rows = r(n_rows)
test_assert "[diagnostics command] target table has one row per target" `diag_rows' == 5
test_assert "[diagnostics command] target scope is recorded" "`diag_scope'" == "targets"
test_assert "[diagnostics command] full target table includes coordinates and diagnostics" colsof(D_post) == 22

quietly cate_dnn_diagnostics, allcandidates
matrix C_post = r(diagnostics)
local candidate_scope "`r(scope)'"
test_assert "[diagnostics command] candidate scope is recorded" "`candidate_scope'" == "candidates"
test_assert "[diagnostics command] paper grid reports four candidates" rowsof(C_post) == 4
test_assert "[diagnostics command] candidate table has documented schema" colsof(C_post) == 13

quietly cate_dnn_diagnostics, candidates(2 4)
test_assert "[diagnostics command] candidate subsetting works" r(n_rows) == 2

tempfile diagnostic_file
quietly cate_dnn_diagnostics, saving("`diagnostic_file'") replace
capture confirm file "`diagnostic_file'"
test_assert "[diagnostics command] saving() writes a Stata dataset" _rc == 0
matrix B_after_diag = e(b)
mata: st_numscalar("_post_diag_b", max(abs(st_matrix("B_post") - st_matrix("B_after_diag"))))
test_assert "[diagnostics command] active estimates are preserved" scalar(_post_diag_b) < 1e-12

test_section "cate_dnn postestimation: stored-result plot builders"

foreach p in estimates diagnostics tuning {
    capture noisily cate_dnn_plot, type(`p') name(cate_post_`p') replace nodraw
    local rc_`p' = _rc
    if !_rc {
        capture graph describe cate_post_`p'
        local graph_rc_`p' = _rc
    }
    else local graph_rc_`p' = 1
    test_assert "[plot] `p' builder runs" `rc_`p'' == 0
    test_assert "[plot] `p' builder creates named graph" `graph_rc_`p'' == 0
}

capture noisily cate_dnn_plot, type(estimates) allcandidates ///
    name(cate_post_est_candidates) replace nodraw
test_assert "[plot] candidate estimate comparison runs" _rc == 0
capture noisily cate_dnn_plot, type(diagnostics) allcandidates ///
    name(cate_post_diag_candidates) replace nodraw
test_assert "[plot] candidate diagnostic comparison runs" _rc == 0

foreach v in mass propensity {
    capture noisily cate_dnn_plot, type(overlap) view(`v') ///
        name(cate_post_overlap_`v') replace nodraw
    test_assert "[plot] overlap `v' view runs" _rc == 0
}

test_section "cate_dnn postestimation: weight and influence builders"

foreach v in profile cumulative kernel {
    capture noisily cate_dnn_plot, type(weights) view(`v') ///
        name(cate_post_weights_`v') replace nodraw
    test_assert "[plot] weight `v' view runs" _rc == 0
}
capture noisily cate_dnn_plot, type(weights) view(profile) allcandidates ///
    name(cate_post_weights_candidates) replace nodraw
test_assert "[plot] candidate weight profiles run" _rc == 0

capture noisily cate_dnn_plot, type(influence) targets(1 3 5) nlabel(2) ///
    name(cate_post_influence) replace nodraw
test_assert "[plot] influence view runs" _rc == 0
if !_rc local influence_n = r(n_plotted)
else local influence_n = .
test_assert "[plot] influence view records selected targets" `influence_n' == 3

capture noisily cate_dnn_plot, type(dashboard) name(cate_post_dashboard) ///
    replace nodraw
test_assert "[plot] dashboard composition runs" _rc == 0
if !_rc local dashboard_n = r(n_plotted)
else local dashboard_n = .
test_assert "[plot] dashboard contains multiple panels" `dashboard_n' >= 4

tempfile plot_stub
local plot_png "`plot_stub'.png"
capture noisily cate_dnn_plot, type(estimates) name(cate_post_export) ///
    replace nodraw saving("`plot_png'")
local export_rc = _rc
capture confirm file "`plot_png'"
local export_file_rc = _rc
test_assert "[plot] graph export runs" `export_rc' == 0
test_assert "[plot] graph export writes PNG" `export_file_rc' == 0

matrix B_after_plots = e(b)
mata: st_numscalar("_post_plot_b", max(abs(st_matrix("B_post") - st_matrix("B_after_plots"))))
test_assert "[plot] active estimates are preserved" scalar(_post_plot_b) < 1e-12
test_assert "[plot] user dataset observation count is preserved" _N == 120

test_section "cate_dnn postestimation: errors"

capture noisily cate_dnn_plot, type(estimates) style(bad)
test_assert "[plot errors] invalid estimate style is rejected" _rc == 198
capture noisily cate_dnn_plot, type(weights) view(kernel) allcandidates
test_assert "[plot errors] candidate kernel comparison is rejected" _rc == 198

preserve
drop x1
capture noisily cate_dnn_plot, type(weights) view(kernel)
local stale_rc = _rc
restore
test_assert "[plot errors] stale observation-level data are rejected" `stale_rc' == 498

quietly cate_dnn y, xcols(x1 x2) x0matrix(X0_post) ///
    pseudooutcome k(2) s(10)
capture noisily cate_dnn_plot, type(overlap) view(mass)
test_assert "[plot errors] overlap rejects precomputed-score fits" _rc == 498
capture noisily cate_dnn_plot, type(tuning)
test_assert "[plot errors] tuning requires a grid fit" _rc == 498
capture noisily cate_dnn_diagnostics, allcandidates
test_assert "[diagnostics errors] candidates require a grid fit" _rc == 498

ereturn clear
capture noisily cate_dnn_plot
test_assert "[plot errors] active cate_dnn estimates are required" _rc == 301
capture noisily cate_dnn_diagnostics
test_assert "[diagnostics errors] active cate_dnn estimates are required" _rc == 301

capture graph drop cate_post_estimates cate_post_diagnostics cate_post_tuning ///
    cate_post_est_candidates cate_post_diag_candidates ///
    cate_post_overlap_mass cate_post_overlap_propensity ///
    cate_post_weights_profile cate_post_weights_cumulative ///
    cate_post_weights_kernel cate_post_weights_candidates ///
    cate_post_influence cate_post_dashboard cate_post_export
