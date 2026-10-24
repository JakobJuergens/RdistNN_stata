version 17
do "tests/test_framework.do"
test_section "fresh and restored diagnostics share matrix-input contracts"
clear
set obs 24
gen double dx=(_n-12)/10
gen double dx2=sin(_n/3)
gen double dy=1+dx^2+.1*sin(_n)
gen byte da=mod(_n,2)
gen double de=.4+.05*dx
gen double dm0=dx
gen double dm1=dx+1
gen double dd=.8+.05*dx
gen double dzero=0
gen double done=1
matrix DC_TARGETS=(-.4,0\ .3,.1\-.4,0)
foreach route in np cate_aipw cate_supplied clate_iv clate_paired {
    foreach method in exact asymp {
        foreach distance in l2 l1 l1_scaled l2_scaled mahalanobis {
            local cmd "np_dnn"
            local extra ""
            if substr("`route'",1,4)=="cate" {
                local cmd "cate_dnn"
                if "`route'"=="cate_aipw" local extra "tcol(da) ecol(de) m0col(dm0) m1col(dm1)"
                else local extra "pseudooutcome"
            }
            if substr("`route'",1,5)=="clate" {
                local cmd "clate_dnn"
                if "`route'"=="clate_iv" local extra "dcol(da) acol(da) ecol(de) my0col(dm0) my1col(dm1) md0col(dzero) md1col(done)"
                else local extra "pseudooutcome denominator(dd)"
            }
            quietly `cmd' dy, xcols(dx dx2) x0matrix(DC_TARGETS) method(`method') distance(`distance') ///
                tuningmode(grid) sgrid(3 5 3 24) sdefault(5) `extra'
            matrix DC_FRESH=e(diagnostics)
            matrix DC_RANKS=e(model_ranks)
            matrix DC_DISTANCE=e(model_distances)
            tempname parent
            estimates store `parent'
            quietly rdistnn_predict, x0matrix(DC_TARGETS) sefit post
            mata: st_numscalar("DC_SAME",all(st_matrix("DC_FRESH"):==st_matrix("e(diagnostics)")) & all(st_matrix("DC_RANKS"):==st_matrix("e(model_ranks)")) & all(st_matrix("DC_DISTANCE"):==st_matrix("e(model_distances)")))
            test_assert "`route' `method' `distance': restored diagnostic values, unavailable masks and frozen ranks are exact" DC_SAME
            quietly estimates restore `parent'
            quietly rdistnn_candidate, candidate(2) post
            mata: st_numscalar("DC_SAME",all(st_matrix("DC_FRESH"):==st_matrix("e(diagnostics)")))
            test_assert "`route' `method' `distance': default candidate retains all diagnostic values and masks" DC_SAME
            estimates drop `parent'
        }
    }
}


test_section "CATE and CLATE retain singleton arm ESS through saved reuse"
clear
set obs 600
gen double ess_x=_n
gen double ess_y=0
gen byte ess_arm=(_n==301)
gen double ess_e=.5
gen double ess_zero=0
gen double ess_one=1
foreach family in cate clate {
    foreach flip in 0 1 {
        replace ess_arm=cond(`flip',_n!=301,_n==301)
        local opts "tcol(ess_arm) ecol(ess_e) m0col(ess_zero) m1col(ess_zero)"
        if "`family'"=="clate" local opts "dcol(ess_arm) acol(ess_arm) ecol(ess_e) my0col(ess_zero) my1col(ess_zero) md0col(ess_zero) md1col(ess_one)"
        quietly `family'_dnn ess_y, xcols(ess_x) x0(0) k(1) method(exact) ///
            tuningmode(grid) sgrid(299 300) sdefault(300) nose `opts'
        local masscol=13+`flip'
        local esscol=22+`flip'
        test_assert "`family' arm `flip': positive rank301 mass has ESS one" ///
            el(e(diagnostics),1,`masscol')>0 & el(e(diagnostics),1,`esscol')==1
        tempname ess_parent
        estimates store `ess_parent'
        quietly rdistnn_candidate, candidate(2) post
        test_assert "`family' arm `flip': default candidate retains ESS one" ///
            el(e(diagnostics),1,`esscol')==1
        quietly estimates restore `ess_parent'
        tempfile ess_saved
        quietly estimates save "`ess_saved'", replace
        preserve
        clear
        mata: mata clear
        quietly estimates use "`ess_saved'"
        quietly rdistnn_predict, x0(0) post
        test_assert "`family' arm `flip': cold restored prediction retains ESS one" ///
            el(e(diagnostics),1,`masscol')>0 & el(e(diagnostics),1,`esscol')==1 & _N==0
        restore
        estimates drop `ess_parent'
    }
}
