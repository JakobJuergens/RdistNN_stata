// test-cate-dnn-tuning.do
// Integrated tuning-grid behavior and deterministic R-package parity.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "cate_dnn: [R parity] paper tuning grid"

clear
input double x1 x2 rho_grid
-1  0   -1.5
 0 -1    0.25
 0  0    1
 1  0    2.5
 0  1   -0.75
 1  1    3
-1 -1    0.5
 2 -0.5  4.25
end
matrix X0_grid = (0.2,-0.1 \ -0.5,0.75 \ 1.2,0.4)

quietly cate_dnn rho_grid, xcols(x1 x2) x0matrix(X0_grid) ///
    k(2) method(exact) pseudooutcome tuningmode(grid)
matrix C_grid = e(tuning_candidates)
matrix E_grid = e(tuning_estimates)
matrix S_grid = e(tuning_se)
matrix D_grid = e(tuning_diagnostics)
matrix B_grid = e(b)
matrix V_grid = e(V)

matrix C_r = (1,1,0,1,0,1,.75,3 \ ///
              2,2,1,1,0,1,1,4 \ ///
              3,3,0,1,0,1,1.25,5 \ ///
              4,4,0,1,0,1,1.5,7)
matrix E_r = ( ///
    .700892857142857, .625, 1.79464285714286 \ ///
    .778571428571428, .364285714285715, 2.38214285714286 \ ///
    1.00892857142857, -.042410714285714, 2.75669642857143 \ ///
    1.5625, -.859375, 2.9375)
matrix S_r = ( ///
    .625263337379882, .776884230206057, .903008364935152 \ ///
    .749642772068327, 1.06502850397805, .885851083905835 \ ///
    .956066158798647, 1.30950711062134, .713367853352273 \ ///
    .992156741649222, 1.24226120064381, .875)
matrix D_r = ( ///
    6.30985915492958, .1875, 5, .848214285714286, 0 \ ///
    4.86593843098312, .25, 5, .964285714285714, 0 \ ///
    3.8478527607362, .3125, 5, 1, 0 \ ///
    2.50980392156863, .4375, 5, 1, 0)

mata {
    st_numscalar("_grid_c_r", max(abs(st_matrix("C_grid") - st_matrix("C_r"))))
    st_numscalar("_grid_e_r", max(abs(st_matrix("E_grid") - st_matrix("E_r"))))
    st_numscalar("_grid_s_r", max(abs(st_matrix("S_grid") - st_matrix("S_r"))))
    st_numscalar("_grid_d_r", max(abs(st_matrix("D_grid") - st_matrix("D_r"))))
    st_numscalar("_grid_default_b", max(abs(st_matrix("B_grid") - st_matrix("E_r")[2,.])))
    st_numscalar("_grid_default_s", max(abs(sqrt(diagonal(st_matrix("V_grid")))' - st_matrix("S_r")[2,.])))
}
test_assert "[R parity] candidate schema values match RdistNN" scalar(_grid_c_r) < 1e-12
test_assert "[R parity] all grid estimates match RdistNN" scalar(_grid_e_r) < 1e-10
test_assert "[R parity] all grid SEs match RdistNN" scalar(_grid_s_r) < 1e-10
test_assert "[R parity] all grid diagnostics match RdistNN" scalar(_grid_d_r) < 1e-10
test_assert "[grid] e(b) is the marked default candidate" scalar(_grid_default_b) < 1e-10
test_assert "[grid] e(V) diagonal is the marked default candidate" scalar(_grid_default_s) < 1e-10
test_assert "[grid] default paper order is s=4" e(s) == 4
test_assert "[grid] default candidate is row 2" e(tuning_default_candidate) == 2
test_assert "[grid] all four candidates are feasible" e(tuning_n_feasible) == 4
test_assert "[grid] all four computed orders are unique" e(tuning_n_unique) == 4
test_assert "[grid] R-aligned tuning mode is recorded" "`e(tuning_mode)'" == "grid"

test_section "cate_dnn: [grid] scalar equivalence and duplicate orders"

forvalues g = 1/4 {
    local sg = C_grid[`g', 8]
    quietly cate_dnn rho_grid, xcols(x1 x2) x0matrix(X0_grid) ///
        k(2) s(`sg') method(exact) pseudooutcome
    matrix B_scalar = e(b)
    matrix V_scalar = e(V)
    mata {
        st_numscalar("_grid_scalar_b", max(abs(st_matrix("B_scalar") - st_matrix("E_r")[`g',.])))
        st_numscalar("_grid_scalar_s", max(abs(sqrt(diagonal(st_matrix("V_scalar")))' - st_matrix("S_r")[`g',.])))
    }
    test_assert "[grid] candidate `g' estimate equals scalar fit" scalar(_grid_scalar_b) < 1e-10
    test_assert "[grid] candidate `g' SE equals scalar fit" scalar(_grid_scalar_s) < 1e-10
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

test_section "cate_dnn: [grid] infeasible rows and validation"

quietly cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) pseudooutcome ///
    tuningmode(grid) sgrid(1 4 8) sdefault(4)
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
    pseudooutcome lambdagrid(.75 1)
test_assert "[errors] lambdagrid() requires grid mode" _rc == 198
capture noisily cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) ///
    pseudooutcome tuningmode(grid) lambdagrid(.75 1) sgrid(3 4)
test_assert "[errors] lambda and s grids are mutually exclusive" _rc == 198
capture noisily cate_dnn rho_grid, xcols(x1 x2) x0(0 0) k(2) ///
    pseudooutcome tuningmode(grid) sgrid(1 3) sdefault(1)
test_assert "[errors] infeasible default candidate fails" _rc == 198
