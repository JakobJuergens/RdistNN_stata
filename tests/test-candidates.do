version 17
do "tests/helpers/fitted-data.do"
test_section "candidate views and prediction preserve the caller"
quietly np_dnn fr_y, xcols(fr_x) x0(0) tuningmode(grid) sgrid(2 3 4) sdefault(3)
tempname parent candidate_estimate expected_b expected_V
estimates store `parent'
matrix FR_PARENT=e(b)
matrix `candidate_estimate'=e(tuning_estimates)
capture noisily rdistnn_candidate, candidate(3)
local candidate_rc=_rc
test_assert "candidate extraction is available" `candidate_rc'==0
if `candidate_rc'==0 {
    matrix FR_CANDIDATE=r(fit)
    mata: st_numscalar("_fr_parent",max(abs(st_matrix("FR_PARENT")-st_matrix("e(b)"))))
    test_assert "candidate extraction preserves active e()" scalar(_fr_parent)==0 & e(s)==3
    test_assert "candidate extraction uses the requested order" abs(FR_CANDIDATE[1,1]-`candidate_estimate'[3,1])<1e-12
    rdistnn_candidate, candidate(3) post
    test_assert "explicit candidate posting creates a reusable fit" e(s)==4 & e(model_version)==3 & "`e(cmd)'"=="np_dnn"
}
matrix FR_TARGETS=(0\.5\2)
quietly np_dnn fr_y, xcols(fr_x) x0matrix(FR_TARGETS) s(3)
matrix `expected_b'=e(b)
matrix `expected_V'=e(V)
quietly estimates restore `parent'
preserve
clear
capture noisily rdistnn_predict, x0matrix(FR_TARGETS) sefit
local prediction_rc=_rc
test_assert "prediction can use a saved fit without the active observations" `prediction_rc'==0
if `prediction_rc'==0 {
    matrix FR_FIT=r(fit)
    matrix FR_V=r(vcov)
    mata: st_numscalar("_fr_prediction",max((max(abs(st_matrix("FR_FIT")-st_matrix("`expected_b'"))),max(abs(st_matrix("FR_V")-st_matrix("`expected_V'"))))))
    test_assert "new-target predictions and covariance match a fresh fit" scalar(_fr_prediction)<1e-11
    test_assert "prediction preserves an empty caller dataset" _N==0
    test_assert "prediction preserves the active order and target batch" e(s)==3 & e(m)==1
}
restore

test_section "standalone candidates across scalar estimators"
foreach family in np cate clate {
    local opts "xcols(fr_x) x0matrix(FR_NEW) tuningmode(grid) sgrid(3 4 4 5) sdefault(4)"
    if "`family'"=="np" quietly np_dnn fr_y, `opts'
    if "`family'"=="cate" quietly cate_dnn fr_y, pseudooutcome `opts'
    if "`family'"=="clate" quietly clate_dnn fr_y, pseudooutcome denominator(fr_den) `opts'
    matrix FR_GRID=e(tuning_estimates)
    quietly rdistnn_candidate, candidate(4) store(fr_candidate)
    matrix FR_FIT=r(fit)
    mata: st_numscalar("_fr_error",max(abs(st_matrix("FR_FIT")-st_matrix("FR_GRID")[4,])))
    test_assert "`family': nondefault candidate matches retained grid estimates" scalar(_fr_error)<1e-12 & e(s)==4
    quietly rdistnn_candidate, candidate(2)
    matrix FR_DUP=r(fit)
    quietly rdistnn_candidate, candidate(3)
    matrix FR_FIT=r(fit)
    mata: st_numscalar("_fr_error",max(abs(st_matrix("FR_DUP")-st_matrix("FR_FIT"))))
    test_assert "`family': duplicate resolved orders retain numerical parity" scalar(_fr_error)==0
    quietly estimates restore fr_candidate
    local fields : e(matrices)
    local gridfield "tuning_candidates"
    local nested : list gridfield in fields
    test_assert "`family': derived candidate is standalone and retains requested identity" e(s)==5 & e(candidate_id)==4 & !`nested'
    quietly `family'_dnn_diagnostics
    quietly `family'_dnn_plot, type(estimates) name(fr_candidate_plot) replace nodraw
    quietly rdistnn_predict, x0(.5) sefit
    test_assert "`family': candidate supports reporting, plots and subsequent prediction" e(model_version)==3 & e(s)==5
    estimates drop fr_candidate
}
