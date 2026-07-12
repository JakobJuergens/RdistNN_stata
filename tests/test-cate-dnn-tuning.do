// Integrated tuning-grid behavior and deterministic R-package parity.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
quietly do "tests/fixtures/r_priority_fixtures.do"

test_section "cate_dnn: [R parity] paper tuning grid"

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
gen double rho_grid = 1 + 0.3*x1 - 0.2*x2 + (-1)^fixture_i*0.25
matrix X0_grid = (0.2,-0.1 \ -0.5,0.75 \ 1.2,0.4)

quietly cate_dnn rho_grid, xcols(x1 x2) x0matrix(X0_grid) ///
    k(2) method(exact) pseudooutcome tuningmode(grid)
matrix C_grid = e(tuning_candidates)
matrix E_grid = e(tuning_estimates)
matrix S_grid = e(tuning_se)
matrix D_grid = e(tuning_diagnostics)
matrix B_grid = e(b)
matrix V_grid = e(V)
matrix ESS_grid = e(ess_targets)

mata {
    _Cg = st_matrix("C_grid"); _Cr = st_matrix("R_cate_grid_candidates")
    _Eg = st_matrix("E_grid"); _Er = st_matrix("R_cate_grid_estimates")
    _Sg = st_matrix("S_grid"); _Sr = st_matrix("R_cate_grid_se")
    _Dg = st_matrix("D_grid"); _Dr = st_matrix("R_cate_grid_diagnostics")
    st_numscalar("_grid_c_r", max(abs(_Cg - _Cr)))
    st_numscalar("_grid_e_r", max(abs(editmissing(_Eg,0) - editmissing(_Er,0))))
    st_numscalar("_grid_e_missing", sum(vec(missing(_Eg) :!= missing(_Er))))
    st_numscalar("_grid_s_r", max(abs(editmissing(_Sg,0) - editmissing(_Sr,0))))
    st_numscalar("_grid_s_missing", sum(vec(missing(_Sg) :!= missing(_Sr))))
    st_numscalar("_grid_d_r", max(abs(editmissing(_Dg,0) - editmissing(_Dr,0))))
    st_numscalar("_grid_d_missing", sum(vec(missing(_Dg) :!= missing(_Dr))))
    st_numscalar("_grid_default_b", max(abs(st_matrix("B_grid") - _Er[2,.])))
    st_numscalar("_grid_default_s", max(abs(sqrt(diagonal(st_matrix("V_grid")))' - _Sr[2,.])))
    st_numscalar("_grid_default_ess", max(abs(st_matrix("ESS_grid") - st_matrix("R_cate_grid_default_ess")')))
}
test_assert "[R parity] candidate order, lambda, and resolved s match RdistNN" scalar(_grid_c_r) < 1e-12
test_assert "[R parity] all available grid estimates match RdistNN" scalar(_grid_e_r) < 1e-11 & scalar(_grid_e_missing) == 0
test_assert "[R parity] all available grid SEs match RdistNN" scalar(_grid_s_r) < 1e-11 & scalar(_grid_s_missing) == 0
test_assert "[R parity] all grid diagnostics and missing patterns match RdistNN" scalar(_grid_d_r) < 1e-11 & scalar(_grid_d_missing) == 0
test_assert "[R parity] default target ESS matches RdistNN" scalar(_grid_default_ess) < 1e-11
test_assert "[grid] e(b) is the marked default candidate" scalar(_grid_default_b) < 1e-11
test_assert "[grid] e(V) diagonal is the marked default candidate" scalar(_grid_default_s) < 1e-11
test_assert "[grid] finite-sample default is lambda=0.4 and s=2" abs(e(lambda)-.4) < 1e-12 & e(s) == 2
test_assert "[grid] default candidate remains row 2" e(tuning_default_candidate) == 2
test_assert "[grid] lambda=0.3 row is retained as infeasible" C_grid[1,4] == 0 & C_grid[1,5] == 1
test_assert "[grid] three candidates are feasible and unique" e(tuning_n_feasible) == 3 & e(tuning_n_unique) == 3
test_assert "[grid] R-aligned tuning mode is recorded" "`e(tuning_mode)'" == "grid"

test_section "cate_dnn: [grid] scalar equivalence and duplicate orders"

forvalues g = 2/4 {
    local sg = C_grid[`g', 8]
    quietly cate_dnn rho_grid, xcols(x1 x2) x0matrix(X0_grid) ///
        k(2) s(`sg') method(exact) pseudooutcome
    matrix B_scalar = e(b)
    matrix V_scalar = e(V)
    mata {
        st_numscalar("_grid_scalar_b", max(abs(st_matrix("B_scalar") - st_matrix("R_cate_grid_estimates")[`g',.])))
        st_numscalar("_grid_scalar_s", max(abs(sqrt(diagonal(st_matrix("V_scalar")))' - st_matrix("R_cate_grid_se")[`g',.])))
    }
    test_assert "[grid] candidate `g' estimate equals scalar fit" scalar(_grid_scalar_b) < 1e-11
    test_assert "[grid] candidate `g' SE equals scalar fit" scalar(_grid_scalar_s) < 1e-11
}

quietly cate_dnn rho_grid, xcols(x1 x2) x0matrix(X0_grid) ///
    k(2) pseudooutcome tuningmode(grid) sgrid(3 4 4 5) sdefault(4)
matrix C_dup = e(tuning_candidates)
matrix E_dup = e(tuning_estimates)
matrix S_dup = e(tuning_se)
matrix D_dup = e(tuning_diagnostics)
mata {
    st_numscalar("_grid_dup_e", max(abs(st_matrix("E_dup")[2,.] - st_matrix("E_dup")[3,.])))
    st_numscalar("_grid_dup_s", max(abs(st_matrix("S_dup")[2,.] - st_matrix("S_dup")[3,.])))
    st_numscalar("_grid_dup_d", max(abs(st_matrix("D_dup")[2,.] - st_matrix("D_dup")[3,.])))
}
test_assert "[grid] duplicate s shares computed_id" C_dup[2,2] == C_dup[3,2]
test_assert "[grid] duplicate s preserves separate candidate rows" rowsof(C_dup) == 4
test_assert "[grid] duplicate estimates are identical" scalar(_grid_dup_e) < 1e-12
test_assert "[grid] duplicate SEs are identical" scalar(_grid_dup_s) < 1e-12
test_assert "[grid] duplicate diagnostics are identical" scalar(_grid_dup_d) < 1e-12
test_assert "[grid] duplicate orders count once in n_unique" e(tuning_n_unique) == 3

quietly cate_dnn rho_grid, xcols(x1 x2) x0(0.2 -0.1) k(1) ///
    pseudooutcome lambda(1)
test_assert "[override] lambda(1) retains the aggressive benchmark" abs(e(lambda)-1) < 1e-12 & e(s) == 6
quietly cate_dnn rho_grid, xcols(x1 x2) x0(0.2 -0.1) k(1) ///
    pseudooutcome tuningmode(grid) lambda(1) lambdagrid(.4 1)
matrix C_lambda_override = e(tuning_candidates)
test_assert "[override] custom lambdagrid() order is preserved" C_lambda_override[1,7] == .4 & C_lambda_override[2,7] == 1
test_assert "[override] lambda(1) selects its custom-grid row" e(tuning_default_candidate) == 2 & e(s) == 6

test_section "cate_dnn: [grid] infeasible rows and validation"

quietly cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) pseudooutcome ///
    tuningmode(grid) sgrid(1 4 12) sdefault(4)
matrix C_skip = e(tuning_candidates)
matrix E_skip = e(tuning_estimates)
test_assert "[grid] k>s row is retained as infeasible" C_skip[1,4] == 0 & C_skip[1,5] == 1
test_assert "[grid] jackknife s=n row is retained as infeasible" C_skip[3,4] == 0 & C_skip[3,5] == 2
test_assert "[grid] infeasible estimate is missing" missing(E_skip[1,1]) & missing(E_skip[3,1])
test_assert "[grid] skipped count includes infeasible rows" e(tuning_n_skipped) == 2
test_assert "[grid] n_unique counts feasible computed orders" e(tuning_n_unique) == 1

capture noisily cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) s(4) ///
    pseudooutcome tuningmode(grid)
test_assert "[errors] tuningmode(grid) rejects s()" _rc == 198
capture noisily cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) ///
    pseudooutcome lambdagrid(.3 .4)
test_assert "[errors] lambdagrid() requires grid mode" _rc == 198
capture noisily cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) ///
    pseudooutcome tuningmode(grid) lambdagrid(.3 .4) sgrid(3 4)
test_assert "[errors] lambda and s grids are mutually exclusive" _rc == 198
capture noisily cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) ///
    pseudooutcome tuningmode(grid) sgrid(1 3) sdefault(1)
test_assert "[errors] infeasible default candidate fails" _rc == 198
