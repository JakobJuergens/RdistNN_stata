// test-mata-jackknife-golden.do
// Golden-value checks for rdistnn_jk_cov_from_ranks().
//
// These mirror the R package's small hand-computed rank-jackknife tests, using
// the Stata rank-shift covariance helper that backs np_dnn/cate_dnn/plqr_dnn.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "rdistnn_jk_cov_from_ranks: [golden] small exact case"

mata {
    _ranks = (1 \ 2 \ 3 \ 4 \ 5)
    _signal = (1 \ 2 \ 3 \ 4 \ 5)
    _V = rdistnn_jk_cov_from_ranks(_ranks, _signal, 1, 2, "exact")
    st_numscalar("_jk_golden_exact", abs(_V[1,1] - 26/45))
}
test_assert "[golden] exact scalar covariance matches hand value" scalar(_jk_golden_exact) < 1e-12

mata {
    _ranks = (1 \ 2 \ 3 \ 4 \ 5)
    _signal = ((1 \ 2 \ 3 \ 4 \ 5), 2 :* (1 \ 2 \ 3 \ 4 \ 5))
    _V = rdistnn_jk_cov_from_ranks(_ranks, _signal, 1, 2, "exact")
    _Vexp = (26/45, 52/45 \ 52/45, 104/45)
    st_numscalar("_jk_golden_vec", max(abs(_V - _Vexp)))
}
test_assert "[golden] exact vector covariance preserves cross-covariances" scalar(_jk_golden_vec) < 1e-12

test_section "rdistnn_jk_cov_from_ranks: [errors] malformed inputs"

capture mata: _V_bad_rows = rdistnn_jk_cov_from_ranks((1 \ 2 \ 3), (1 \ 2), 1, 2, "exact")
test_assert "[errors] signal/rank row mismatch rejected" _rc != 0

capture mata: _V_bad_method = rdistnn_jk_cov_from_ranks((1 \ 2 \ 3), (1 \ 2 \ 3), 1, 2, "bad")
test_assert "[errors] unknown jackknife method rejected" _rc != 0

test_section "rdistnn_jk_cov_from_ranks: [det] deterministic repeated calls"

mata {
    rseed(721)
    _X = rnormal(24, 2, 0, 1)
    _ranks = rdistnn_ranks(_X, (0.2, -0.3), "l2")
    _signal = rnormal(24, 3, 0, 1)
    _V1 = rdistnn_jk_cov_from_ranks(_ranks, _signal, 2, 5, "asymp")
    _V2 = rdistnn_jk_cov_from_ranks(_ranks, _signal, 2, 5, "asymp")
    st_numscalar("_jk_repeat_diff", max(abs(_V1 - _V2)))
}
test_assert "[det] repeated asymp vector covariance calls match exactly" scalar(_jk_repeat_diff) == 0
