// test-weights-dnn-smoke.do
// Tests for the weights_dnn command: basic behavior, return values, and metadata.
//
// Tags: [smoke] [meta] [edge]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_basic, n(200) seed(42)

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [smoke] basic run"
// ════════════════════════════════════════════════════════════

// [smoke] command runs without error and returns expected scalars
weights_dnn x1 x2, x0(0 0) k(5) generate(w_test)
local neff = r(neff)
local n_r  = r(n)
local k_r  = r(k)
local s_r  = r(s)
local kx_r = r(kx)

test_assert "[smoke] r(neff) > 1"              `neff' > 1
test_assert "[smoke] r(neff) <= n"             `neff' <= `n_r'
test_assert "[smoke] r(n) = 200"               `n_r'  == 200
test_assert "[smoke] r(k) = 5"                 `k_r'  == 5
test_assert "[smoke] r(kx) = 2"                `kx_r' == 2
test_assert "[smoke] r(s) > 0"                 `s_r'  > 0

// [smoke] generated weights sum to 1
quietly summarize w_test
local wsum = r(sum)
test_assert "[smoke] generate: weights sum to 1" abs(`wsum' - 1) < 1e-8

// [smoke] all generated weights are non-negative
quietly summarize w_test
local wmin = r(min)
test_assert "[smoke] generate: weights non-negative" `wmin' >= 0

// [smoke] generated variable covers all observations in sample
quietly count if w_test == .
local nmiss = r(N)
test_assert "[smoke] generate: no missing values in sample" `nmiss' == 0

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [meta] returned scalars"
// ════════════════════════════════════════════════════════════

// [meta] explicit k and s are passed through correctly
weights_dnn x1 x2, x0(0 0) k(3) s(10)
local k_ret = r(k)
local s_ret = r(s)
test_assert "[meta] explicit k passed through" `k_ret' == 3
test_assert "[meta] explicit s passed through" `s_ret' == 10

// [meta] auto-s default: floor(0.4*n^((kx+1)/(kx+2))) with kx=2, n=200
//   = floor(0.4*200^(3/4)) = floor(21.27...) = 21
weights_dnn x1 x2, x0(0 0) k(3)
local s_auto = r(s)
test_assert "[meta] auto s = floor(0.4*n^(3/4))" `s_auto' == 21

// [meta] k defaults to the paper/R value 1 when omitted.
weights_dnn x1 x2, x0(0 0) s(10)
test_assert "[meta] omitted k defaults to 1" r(k) == 1

// [edge] Weight construction permits the mathematical endpoint s=n.
clear
set obs 12
gen double x1 = _n
weights_dnn x1, x0(0) s(12) method(exact) generate(w_sn)
quietly count if w_sn > 0
local npos_sn = r(N)
quietly summarize w_sn, meanonly
test_assert "[edge] s=n with default k=1 selects only the nearest rank" `npos_sn' == 1
test_assert "[edge] exact s=n weights remain normalized" abs(r(sum) - 1) < 1e-12

drop w_sn
weights_dnn x1, x0(0) k(5) s(5) method(exact) generate(w_ks)
quietly summarize w_ks, meanonly
test_assert "[edge] k=s gives uniform weights over the sample" ///
    abs(r(min) - 1/12) < 1e-12 & abs(r(max) - 1/12) < 1e-12

drop w_ks
weights_dnn x1, x0(0) k(12) s(12) method(exact) generate(w_kn)
quietly summarize w_kn, meanonly
test_assert "[edge] k=s=n gives normalized uniform weights" ///
    abs(r(sum) - 1) < 1e-12 & abs(r(min) - 1/12) < 1e-12

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [edge] single covariate"
// ════════════════════════════════════════════════════════════

// [edge] kx = 1 (one covariate): should run fine without mahal restriction
weights_dnn x1, x0(0) k(3) s(3)
local neff1 = r(neff)
test_assert "[edge] kx=1 runs without error" `neff1' > 0

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [edge] if/in restrictions"
// ════════════════════════════════════════════════════════════

// [edge] if restriction reduces sample correctly
rdnn_dgp_basic, n(200) seed(42)
quietly count if x1 > 0
local n_pos = r(N)

weights_dnn x1 x2 if x1 > 0, x0(0.5 0) k(3)
local n_ret = r(n)
test_assert "[edge] if restriction: r(n) matches subsample size" `n_ret' == `n_pos'
