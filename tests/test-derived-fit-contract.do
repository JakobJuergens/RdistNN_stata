version 17
do "tests/test_framework.do"
test_section "derived fit coefficient stripes and verified estimation samples"
clear
set obs 12
gen double fx=(_n-6)/5
gen double fy=1+fx^2+.2*sin(_n)
gen double fd=.8+.05*fx
matrix FT=(-.4\ .4)
foreach targetset in scalar batch {
    local targetopt "x0(0)"
    local names_expected "target1"
    local contrast "target1"
    if "`targetset'"=="batch" {
        local targetopt "x0matrix(FT)"
        local names_expected "target1 target2"
        local contrast "target1-target2"
    }
    quietly np_tdnn fy, xcols(fx) `targetopt' s1(2) s2(3)
    matrix TD_base_b=e(b)
    matrix TD_base_V=e(V)
    quietly rdistnn_predict, `targetopt' sefit post
    capture quietly lincom `contrast'
    local rc=_rc
    test_assert "TDNN `targetset': prediction preserves existing coefficient contrasts" `rc'==0
    foreach field in b V estimates se_targets estimate_vcov {
        local names : colnames e(`field')
        test_assert "TDNN `targetset': `field' retains target coefficient names" "`names'"=="`names_expected'"
    }
    mata: st_numscalar("td_b_diff",max(abs(st_matrix("e(b)")-st_matrix("TD_base_b")))); st_numscalar("td_V_diff",max(abs(st_matrix("e(V)")-st_matrix("TD_base_V"))))
    test_assert "TDNN `targetset': prediction retains estimates" td_b_diff<1e-12
    test_assert "TDNN `targetset': prediction retains covariance" td_V_diff<1e-12
}
matrix drop TD_base_b TD_base_V
scalar drop td_b_diff td_V_diff
foreach family in np cate clate {
    local extra ""
    local coef "mu"
    if "`family'"=="cate" {
        local extra "pseudooutcome"
        local coef "tau"
    }
    if "`family'"=="clate" {
        local extra "pseudooutcome denominator(fd)"
        local coef "theta"
    }
    quietly `family'_dnn fy, xcols(fx) x0matrix(FT) tuningmode(grid) sgrid(2 3) sdefault(3) `extra'
    quietly rdistnn_candidate, candidate(1) post
    capture quietly lincom `coef'1-`coef'2
    local rc=_rc
    test_assert "`family': candidate supports documented coefficient contrasts" `rc'==0
    quietly rdistnn_predict, x0matrix(FT) sefit post
    capture quietly lincom `coef'1-`coef'2
    local rc=_rc
    test_assert "`family': prediction supports documented coefficient contrasts" `rc'==0
    quietly rdistnn_predict, x0(0) post
    local names : colnames e(b)
    test_assert "`family': scalar prediction retains family coefficient" "`names'"=="`coef'"
}
foreach family in np tdnn cate clate {
    foreach state in unchanged outside inside reordered cleared {
        clear
        set obs 12
        gen double fx=(_n-6)/5
        gen double fy=1+fx^2+.2*sin(_n)
        gen double fd=.8+.05*fx
        gen long fid=_n
        replace fy=. in 5
        local extra "s(3)"
        local cmd "`family'_dnn"
        if "`family'"=="tdnn" {
            local cmd "np_tdnn"
            local extra "s1(2) s2(3)"
        }
        if "`family'"=="cate" local extra "`extra' pseudooutcome"
        if "`family'"=="clate" local extra "`extra' pseudooutcome denominator(fd)"
        quietly `cmd' fy if fid>=3 & fid<=10, xcols(fx) x0matrix(FT) `extra'
        local n=e(N)
        if "`state'"=="outside" quietly replace fy=999 in 1
        if "`state'"=="inside" quietly replace fy=999 in 4
        if "`state'"=="reordered" quietly gsort -fid
        if "`state'"=="cleared" clear
        quietly rdistnn_predict, x0matrix(FT) sefit post
        quietly count if e(sample)
        local expected=cond(inlist("`state'","unchanged","outside"),`n',0)
        test_assert "`family' `state': e(sample) uses verified original rows only" r(N)==`expected'
    }
}
clear
set obs 4
gen double fx=_n
gen double fy=1
gen double fd=(_n>2)
matrix PART=(0\5)
quietly clate_dnn fy, denominator(fd) pseudooutcome xcols(fx) x0matrix(PART) tuningmode(grid) sgrid(3 4) sdefault(3) nose
quietly rdistnn_candidate, candidate(2) post
local names : colnames e(b)
test_assert "partially undefined candidate retains original theta2 stripe" "`names'"=="theta2"
quietly rdistnn_predict, x0matrix(PART) post
local names : colnames e(b)
test_assert "partially undefined prediction retains original theta2 stripe" "`names'"=="theta2"
