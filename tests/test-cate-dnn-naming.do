// test-cate-dnn-naming.do
// Canonical option and result metadata.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "cate_dnn: [naming] canonical AIPW metadata"

rdnn_dgp_basic, n(100) seed(2301)
quietly cate_dnn y, tcol(D) xcols(x1 x2) x0(0.1 -0.2) k(2) s(7) ///
    ecol(e) m0col(m0) m1col(m1)
matrix B_r_names = e(b)
matrix V_r_names = e(V)
matrix ESS_r_names = e(ess_targets)
matrix NEFF_r_names = e(neff_targets)
local x_cols_name "`e(x_cols)'"
local t_col_name "`e(t_col)'"
local e_col_name "`e(e_col)'"
local m0_col_name "`e(m0_col)'"
local m1_col_name "`e(m1_col)'"
local jackknife_name "`e(jackknife)'"
local p_name = e(p)
local ess_name = e(ess)
local neff_name = e(neff)
local ess_min_name = e(ess_min)
local neff_min_name = e(neff_min)
local ess_max_name = e(ess_max)
local neff_max_name = e(neff_max)

mata {
    st_numscalar("_names_ess", max(abs(st_matrix("ESS_r_names") - st_matrix("NEFF_r_names"))))
}
test_assert "[naming] e(ess_targets) aliases e(neff_targets)" scalar(_names_ess) < 1e-12
test_assert "[naming] e(ess) aliases e(neff)" abs(`ess_name' - `neff_name') < 1e-12
test_assert "[naming] e(ess_min) aliases e(neff_min)" abs(`ess_min_name' - `neff_min_name') < 1e-12
test_assert "[naming] e(ess_max) aliases e(neff_max)" abs(`ess_max_name' - `neff_max_name') < 1e-12
test_assert "[naming] e(x_cols) records canonical covariates" "`x_cols_name'" == "x1 x2"
test_assert "[naming] e(t_col) records treatment" "`t_col_name'" == "D"
test_assert "[naming] e(e_col) records propensity" "`e_col_name'" == "e"
test_assert "[naming] e(m0_col) records control nuisance" "`m0_col_name'" == "m0"
test_assert "[naming] e(m1_col) records treated nuisance" "`m1_col_name'" == "m1"
test_assert "[naming] e(jackknife) matches R fast route" "`jackknife_name'" == "fast"
test_assert "[naming] e(p) records covariate dimension" `p_name' == 2
