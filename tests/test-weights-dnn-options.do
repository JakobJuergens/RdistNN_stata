// test-weights-dnn-options.do
// Tests for weights_dnn option combinations: method, distance, generate.
//
// Tags: [opts] [agree] [distance]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_basic, n(300) seed(55)

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [opts] method=asymp"
// ════════════════════════════════════════════════════════════

// [opts] asymp method runs without error
weights_dnn x1 x2, x0(0 0) k(4) s(12) method(asymp) generate(w_ap)
local neff_ap = r(neff)
test_assert "[opts] method=asymp runs and returns neff > 1" `neff_ap' > 1

quietly summarize w_ap
local wsum_ap = r(sum)
test_assert "[opts] method=asymp: weights sum to 1" abs(`wsum_ap' - 1) < 1e-8

// [agree] exact and asymp neff are close for large n (within 5%)
weights_dnn x1 x2, x0(0 0) k(4) s(12) method(exact)
local neff_ex = r(neff)
local neff_ratio = abs(`neff_ap' - `neff_ex') / `neff_ex'
test_assert "[agree] exact vs asymp neff within 5%" `neff_ratio' < 0.05

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [distance] mahalanobis"
// ════════════════════════════════════════════════════════════

// [distance] mahalanobis runs without error
weights_dnn x1 x2, x0(0 0) k(4) s(12) distance(mahalanobis) generate(w_mh)
local neff_mh = r(neff)
test_assert "[distance] mahalanobis runs and neff > 1" `neff_mh' > 1

quietly summarize w_mh
local wsum_mh = r(sum)
test_assert "[distance] mahalanobis: weights sum to 1" abs(`wsum_mh' - 1) < 1e-8

// [distance] mahalanobis weights differ from l2 weights (different neighborhoods)
//   For correlated covariates, the two metrics give different orderings.
weights_dnn x1 x2, x0(0 0) k(4) s(12) distance(l2) generate(w_l2)

// compute max absolute difference between l2 and mahal weights
mata {
    _wl2 = st_data(., "w_l2")
    _wmh = st_data(., "w_mh")
    st_numscalar("_r", max(abs(_wl2 - _wmh)))
}
test_assert "[distance] mahal != l2 weights (different neighborhoods)" scalar(_r) > 0

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [opts] generate without generate (no variable created)"
// ════════════════════════════════════════════════════════════

// [opts] calling without generate() should not create a variable
cap drop w_novar
weights_dnn x1 x2, x0(0 0) k(3)
cap confirm variable w_novar
test_assert "[opts] no generate() => no variable created" _rc != 0

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [agree] generate vs no-generate: neff matches"
// ════════════════════════════════════════════════════════════

// [agree] neff from generate run equals neff from no-generate run
weights_dnn x1 x2, x0(1 0) k(5) s(15)
local neff1 = r(neff)

weights_dnn x1 x2, x0(1 0) k(5) s(15) generate(w_check)
local neff2 = r(neff)
gen double w_check_sq = w_check^2
quietly summarize w_check_sq, meanonly
local neff_from_w = 1 / r(sum)

local diff12 = abs(`neff1' - `neff2')
local diff1w = abs(`neff1' - `neff_from_w')

test_assert "[agree] neff same with and without generate()" `diff12' < 1e-8
test_assert "[agree] r(neff) = 1/sum(w^2)"                  `diff1w' < 1e-6


test_section "weights notice facts use target-specific exact distance groups"
clear
set obs 4
gen double notice_x=_n
matrix NOTICE_TARGETS=(0\2.5)
quietly weights_dnn notice_x, x0matrix(NOTICE_TARGETS) s(1)
matrix NOTICE_FACTS=r(target_notices)
test_assert "quiet weights count tied targets and uniform concentration" r(n_tied_targets)==1 & r(n_weak_targets)==2
test_assert "weights facts distinguish symmetric-distance ties from untied profiles" NOTICE_FACTS[1,1]==0 & NOTICE_FACTS[2,1]==1 & NOTICE_FACTS[1,2]==1 & NOTICE_FACTS[2,2]==1
quietly weights_tdnn notice_x, x0matrix(NOTICE_TARGETS) s1(1) s2(2)
matrix NOTICE_FACTS=r(target_notices)
test_assert "quiet signed weights count only the actual tied target" r(n_tied_targets)==1 & r(n_weak_targets)==0 & NOTICE_FACTS[1,1]==0 & NOTICE_FACTS[2,1]==1
