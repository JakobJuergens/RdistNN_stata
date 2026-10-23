// test-cate-dnn-batch.do
// Batch-target parity, joint covariance, score-path, and replay tests.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "cate_dnn: [batch] x0matrix() matches separate x0() calls"

rdnn_dgp_basic, n(120) seed(2201)
matrix X0_batch = (0,0 \ 0.4,-0.2 \ -0.3,0.5)
cate_dnn y, tcol(D) xcols(x1 x2) x0matrix(X0_batch) k(2) s(8) ///
    ecol(e) m0col(m0) m1col(m1)
matrix B_batch = e(b)
matrix V_batch = e(V)
matrix T_batch = e(targets)
matrix N_batch = e(neff_targets)
local m_batch = e(m)
local source_batch "`e(target_source)'"
local x0matrix_batch "`e(x0matrix)'"

matrix B_single = J(1, 3, .)
matrix S_single = J(1, 3, .)
matrix N_single = J(1, 3, .)
forvalues j = 1/3 {
    local tx1 = X0_batch[`j', 1]
    local tx2 = X0_batch[`j', 2]
    quietly cate_dnn y, tcol(D) xcols(x1 x2) x0(`tx1' `tx2') k(2) s(8) ///
        ecol(e) m0col(m0) m1col(m1)
    matrix B_single[1, `j'] = _b[tau]
    matrix S_single[1, `j'] = _se[tau]
    matrix N_single[1, `j'] = e(neff)
}

mata {
    _B_batch = st_matrix("B_batch")
    _V_batch = st_matrix("V_batch")
    _S_batch = sqrt(diagonal(_V_batch))'
    st_numscalar("_batch_bdiff", max(abs(_B_batch - st_matrix("B_single"))))
    st_numscalar("_batch_sdiff", max(abs(_S_batch - st_matrix("S_single"))))
    st_numscalar("_batch_ndiff", max(abs(st_matrix("N_batch") - st_matrix("N_single"))))
    st_numscalar("_batch_tdiff", max(abs(st_matrix("T_batch") - st_matrix("X0_batch"))))
    st_numscalar("_batch_vsym", max(abs(_V_batch - _V_batch')))
    st_numscalar("_batch_vmineig", min(symeigenvalues(_V_batch)))
}
test_assert "[batch] estimates match separate x0() calls" scalar(_batch_bdiff) < 1e-10
test_assert "[batch] diagonal SEs match separate x0() calls" scalar(_batch_sdiff) < 1e-10
test_assert "[batch] effective sample sizes match separate calls" scalar(_batch_ndiff) < 1e-10
test_assert "[batch] e(targets) preserves target values" scalar(_batch_tdiff) < 1e-12
test_assert "[batch] joint e(V) is symmetric" scalar(_batch_vsym) < 1e-10
test_assert "[batch] joint e(V) is positive semidefinite" scalar(_batch_vmineig) > -1e-8
test_assert "[batch] e(m) records three targets" `m_batch' == 3
test_assert "[batch] target source records X0" "`source_batch'" == "X0"
test_assert "[batch] e(x0matrix) records source matrix" "`x0matrix_batch'" == "X0_batch"

test_section "cate_dnn: [subset identities] deterministic batch exact-distance ties"

clear
input double x1 x2 rho_r_golden
-1  0   -1.5
 0 -1    0.25
 0  0    1
 1  0    2.5
 0  1   -0.75
 1  1    3
-1 -1    0.5
 2 -0.5  4.25
end
matrix X0_r_golden = (0.2,-0.1 \ -0.5,0.75 \ 1.2,0.4)
quietly cate_dnn rho_r_golden, xcols(x1 x2) x0matrix(X0_r_golden) ///
    k(2) s(4) method(exact) pseudooutcome
matrix B_r_golden = e(b)
matrix V_r_golden = e(V)
matrix N_r_golden = e(neff_targets)
do "tests/helpers/distance-reference.do"
mata {
    _golden_D=st_matrix("e(model_distances)")
    _golden_S=st_data(.,"rho_r_golden")
    _golden_W=rdistnn_test_subset_weights(_golden_D,2,4)
    _golden_L=rdistnn_test_subset_deletions(_golden_D,_golden_S,2,4)
    _golden_V=(7/8)*cross(_golden_L:-mean(_golden_L),_golden_L:-mean(_golden_L))
    _tau_r_ref=_golden_S'*_golden_W
    _se_r_ref=sqrt(diagonal(_golden_V))'
    _neff_r_ref=1:/colsum(_golden_W:^2)
    _se_stata = sqrt(diagonal(st_matrix("V_r_golden")))'
    st_numscalar("_batch_r_tau", max(abs(st_matrix("B_r_golden") - _tau_r_ref)))
    st_numscalar("_batch_r_se", max(abs(_se_stata - _se_r_ref)))
    st_numscalar("_batch_r_neff", max(abs(st_matrix("N_r_golden") - _neff_r_ref)))
}
test_assert "[subset identities] batch estimates match independent tied subset selection" scalar(_batch_r_tau) < 1e-10
test_assert "[subset identities] batch marginal SEs match independent subset deletions" scalar(_batch_r_se) < 1e-10
test_assert "[subset identities] batch effective sample sizes use averaged actual weights" scalar(_batch_r_neff) < 1e-10

test_section "cate_dnn: [batch] one-row x0matrix() preserves scalar contract"

rdnn_dgp_basic, n(120) seed(2201)
matrix X0_one = (0.25,-0.15)
quietly cate_dnn y, tcol(D) xcols(x1 x2) x0matrix(X0_one) k(2) s(8) ///
    ecol(e) m0col(m0) m1col(m1)
local tau_x0mat_one = _b[tau]
local se_x0mat_one = _se[tau]
local neff_x0mat_one = e(neff)
local m_x0mat_one = e(m)

quietly cate_dnn y, tcol(D) xcols(x1 x2) x0(0.25 -0.15) k(2) s(8) ///
    ecol(e) m0col(m0) m1col(m1)
test_assert "[batch] one-row x0matrix() estimate matches x0()" abs(`tau_x0mat_one' - _b[tau]) < 1e-10
test_assert "[batch] one-row x0matrix() SE matches x0()" abs(`se_x0mat_one' - _se[tau]) < 1e-10
test_assert "[batch] one-row x0matrix() neff matches x0()" abs(`neff_x0mat_one' - e(neff)) < 1e-10
test_assert "[batch] one-row x0matrix() records e(m)=1" `m_x0mat_one' == 1

test_section "cate_dnn: [batch] Mahalanobis uses common fixed metric"

matrix X0_mahal = (0.2,-0.1 \ -0.4,0.35)
quietly cate_dnn y, tcol(D) xcols(x1 x2) x0matrix(X0_mahal) k(2) s(8) ///
    ecol(e) m0col(m0) m1col(m1) distance(mahalanobis)
matrix B_mahal_batch = e(b)
matrix V_mahal_batch = e(V)
matrix B_mahal_single = J(1, 2, .)
matrix S_mahal_single = J(1, 2, .)
forvalues j = 1/2 {
    local tx1 = X0_mahal[`j', 1]
    local tx2 = X0_mahal[`j', 2]
    quietly cate_dnn y, tcol(D) xcols(x1 x2) x0(`tx1' `tx2') k(2) s(8) ///
        ecol(e) m0col(m0) m1col(m1) distance(mahalanobis)
    matrix B_mahal_single[1, `j'] = _b[tau]
    matrix S_mahal_single[1, `j'] = _se[tau]
}
mata {
    _S_mahal_batch = sqrt(diagonal(st_matrix("V_mahal_batch")))'
    st_numscalar("_batch_mahal_b", max(abs(st_matrix("B_mahal_batch") - st_matrix("B_mahal_single"))))
    st_numscalar("_batch_mahal_s", max(abs(_S_mahal_batch - st_matrix("S_mahal_single"))))
}
test_assert "[batch] Mahalanobis estimates match separate calls" scalar(_batch_mahal_b) < 1e-10
test_assert "[batch] Mahalanobis SEs match separate calls" scalar(_batch_mahal_s) < 1e-10

test_section "cate_dnn: [batch] score paths preserve joint results"

rdnn_dgp_basic, n(100) seed(2202)
matrix X0_scores = (0.1,-0.2 \ -0.4,0.3)
gen double eclip_b = min(max(e, 1e-8), 1 - 1e-8)
gen double rho_aipw_b = D*(y-m1)/eclip_b - (1-D)*(y-m0)/(1-eclip_b) + m1-m0

quietly cate_dnn y, tcol(D) xcols(x1 x2) x0matrix(X0_scores) k(2) s(7) ///
    ecol(e) m0col(m0) m1col(m1)
matrix B_aipw_b = e(b)
matrix V_aipw_b = e(V)
quietly cate_dnn rho_aipw_b, xcols(x1 x2) x0matrix(X0_scores) k(2) s(7) pseudooutcome
matrix B_pseudo_b = e(b)
matrix V_pseudo_b = e(V)

mata {
    st_numscalar("_batch_aipw_b", max(abs(st_matrix("B_aipw_b") - st_matrix("B_pseudo_b"))))
    st_numscalar("_batch_aipw_v", max(abs(st_matrix("V_aipw_b") - st_matrix("V_pseudo_b"))))
}
test_assert "[batch] AIPW and precomputed estimates agree" scalar(_batch_aipw_b) < 1e-10
test_assert "[batch] AIPW and precomputed covariance agree" scalar(_batch_aipw_v) < 1e-10

test_section "cate_dnn: [batch] nose and replay"

matrix X0_nose = (0,0 \ 0.5,0.5)
quietly cate_dnn rho_aipw_b, xcols(x1 x2) x0matrix(X0_nose) s(100) pseudooutcome nose
matrix B_nose_batch = e(b)
matrix N_nose_batch = e(neff_targets)
local se_method_batch "`e(se_method)'"
local m_nose_batch = e(m)
capture matrix list e(V)
local rc_batch_v = _rc

quietly cate_dnn
matrix B_replay_batch = e(b)
mata: st_numscalar("_batch_replay", max(abs(st_matrix("B_nose_batch") - st_matrix("B_replay_batch"))))
local bcols = colsof(B_nose_batch)
local ncols = colsof(N_nose_batch)
test_assert "[batch] nose posts one coefficient per target" `bcols' == 2
test_assert "[batch] nose posts one neff per target" `ncols' == 2
test_assert "[batch] nose records no SE method" "`se_method_batch'" == "none"
test_assert "[batch] nose does not post e(V)" `rc_batch_v' != 0
test_assert "[batch] nose records two targets" `m_nose_batch' == 2
test_assert "[batch] replay preserves coefficients" scalar(_batch_replay) < 1e-12
