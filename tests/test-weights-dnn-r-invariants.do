// test-weights-dnn-r-invariants.do
// Command-level weights_dnn invariants ported from the R package tests.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "weights_dnn: [identity] command weights equal Mata rank helpers"

rdnn_dgp_basic, n(90) seed(801)
weights_dnn x1 x2, x0(0.2 -0.3) k(3) s(8) method(exact) generate(w_cmd_ex)
mata {
    _X = st_data(., "x1 x2")
    _w_cmd = st_data(., "w_cmd_ex")
    _ranks = rdistnn_localization_ranks(rdistnn_distance_matrix(_X,(0.2, -0.3),"l2"))
    _w_ref = rdistnn_weights_from_ranks(_ranks, 3, 8, "exact")
    st_numscalar("_weights_exact_diff", max(abs(_w_cmd - _w_ref)))
}
test_assert "[identity] exact generated weights match Mata reference" scalar(_weights_exact_diff) < 1e-12

weights_dnn x1 x2, x0(0.2 -0.3) k(3) s(8) method(asymp) generate(w_cmd_ap)
mata {
    _X = st_data(., "x1 x2")
    _w_cmd = st_data(., "w_cmd_ap")
    _ranks = rdistnn_localization_ranks(rdistnn_distance_matrix(_X,(0.2, -0.3),"l2"))
    _w_ref = rdistnn_weights_from_ranks(_ranks, 3, 8, "asymp")
    st_numscalar("_weights_asymp_diff", max(abs(_w_cmd - _w_ref)))
}
test_assert "[identity] asymp generated weights match Mata reference" scalar(_weights_asymp_diff) < 1e-10

test_section "weights_dnn: [ties] deterministic plateau weights"

clear
set obs 10
gen double x1 = 0
gen double x2 = 0
weights_dnn x1 x2, x0(0 0) k(2) s(4) method(exact) generate(w_tie)
mata {
    _w_cmd = st_data(., "w_tie")
    _w_ref = J(10,1,1/10)
    st_numscalar("_weights_tie_diff", max(abs(_w_cmd - _w_ref)))
    st_numscalar("_weights_tie_sum", abs(sum(_w_cmd) - 1))
}
test_assert "[ties] all-equal X averages all ranks to sample-mean weights" scalar(_weights_tie_diff) < 1e-12
test_assert "[ties] plateau weights sum to one" scalar(_weights_tie_sum) < 1e-12

test_section "weights_dnn: [ifin] generated weights respect estimation sample"

rdnn_dgp_basic, n(100) seed(802)
gen byte keep_w = (_n >= 5 & _n <= 85 & x1 > -0.25)
quietly count if keep_w
local n_keep = r(N)

weights_dnn x1 x2 if keep_w, x0(0.1 -0.2) k(2) s(6) generate(w_ifin)
local n_ret = r(n)
quietly count if keep_w & w_ifin < .
local n_nonmiss_in = r(N)
quietly count if !keep_w & w_ifin < .
local n_nonmiss_out = r(N)
quietly summarize w_ifin if keep_w, meanonly
local w_ifin_sum = r(sum)

test_assert "[ifin] r(n) matches marked sample" `n_ret' == `n_keep'
test_assert "[ifin] all in-sample weights stored" `n_nonmiss_in' == `n_keep'
test_assert "[ifin] out-of-sample weights remain missing" `n_nonmiss_out' == 0
test_assert "[ifin] restricted weights sum to one" abs(`w_ifin_sum' - 1) < 1e-10

test_section "weights_dnn: [edge] one-dimensional signed target"

rdnn_dgp_basic, n(70) seed(803)
weights_dnn x1, x0(-0.4) k(2) s(7) generate(w_1d)
quietly summarize w_1d, meanonly
local w_1d_sum = r(sum)
local w_1d_min = r(min)
local neff_1d = r(neff)
test_assert "[edge] one-dimensional signed x0 weights sum to one" abs(`w_1d_sum' - 1) < 1e-10
test_assert "[edge] one-dimensional signed x0 weights non-negative" `w_1d_min' >= 0
test_assert "[edge] one-dimensional signed x0 neff returned" `neff_1d' > 1

test_section "weights_dnn: [errors] generate() protects existing variables"

rdnn_dgp_basic, n(50) seed(804)
gen double w_exists = .
capture weights_dnn x1 x2, x0(0 0) k(2) s(6) generate(w_exists)
test_assert "[errors] generate() refuses to overwrite existing variable" _rc != 0
