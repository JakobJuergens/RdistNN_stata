// test-cate-dnn-diagnostics.do
// Per-target weight, distance, and overlap diagnostics with R parity.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "cate_dnn: [R parity] target weight and distance diagnostics"

clear
input double x1 x2 rho_diag
-1  0   -1.5
 0 -1    0.25
 0  0    1
 1  0    2.5
 0  1   -0.75
 1  1    3
-1 -1    0.5
 2 -0.5  4.25
end
matrix X0_diag = (0.2,-0.1 \ -0.5,0.75 \ 1.2,0.4)
quietly cate_dnn rho_diag, xcols(x1 x2) x0matrix(X0_diag) ///
    k(2) s(4) method(exact) pseudooutcome
matrix D_diag = e(diagnostics)
matrix T_diag = e(target_diagnostics)

quietly do "tests/fixtures/r_priority_fixtures.do"
matrix D_r_base = R_cate_diagnostics_base
mata {
    _D_diag = st_matrix("D_diag")
    st_numscalar("_diag_r_base", max(abs(_D_diag[,1..10] - st_matrix("D_r_base")[,1..10])))
    st_numscalar("_diag_alias", max(abs(_D_diag - st_matrix("T_diag"))))
}
test_assert "[R parity] weight and mean-distance diagnostics match RdistNN" scalar(_diag_r_base) < 1e-10
test_assert "[diagnostics] R-style and explicit target aliases match" scalar(_diag_alias) < 1e-12
test_assert "[diagnostics] pseudo-outcome records overlap unavailable" D_diag[1,12] == 0
test_assert "[diagnostics] pseudo-outcome records propensity unavailable" D_diag[1,16] == 0
test_assert "[diagnostics] unavailable overlap values are missing" missing(D_diag[1,13]) & missing(D_diag[1,17])

preserve
replace x1 = .a in 1
quietly cate_dnn rho_diag, xcols(x1 x2) x0(0.2 -0.1) ///
    k(2) s(4) method(exact) pseudooutcome
matrix D_complete_retained = e(diagnostics)
local n_complete_retained = e(N)
restore
test_assert "[diagnostics] incomplete raw covariate is excluded before retention" ///
    `n_complete_retained' == 7
test_assert "[diagnostics] complete retained sample yields finite distance radii" ///
    D_complete_retained[1,11] < . & D_complete_retained[1,20] < .

test_section "cate_dnn: [R parity] treatment overlap diagnostics"

gen double y = cond(_n==1,-.5,cond(_n==2,1.2,cond(_n==3,.4, ///
    cond(_n==4,2.1,cond(_n==5,.7,cond(_n==6,1.8,cond(_n==7,-.2,2.7)))))))
gen byte D = inlist(_n,2,4,6,7)
gen double e = cond(_n==1,.2,cond(_n==2,.35,cond(_n==3,.4, ///
    cond(_n==4,.55,cond(_n==5,.6,cond(_n==6,.7,cond(_n==7,.3,.8)))))))
gen double m0 = cond(_n==1,-.2,cond(_n==2,.1,cond(_n==3,.3, ///
    cond(_n==4,.6,cond(_n==5,.4,cond(_n==6,.7,cond(_n==7,-.1,1)))))))
gen double m1 = m0 + 1

quietly cate_dnn y, xcols(x1 x2) x0matrix(X0_diag) tcol(D) ecol(e) ///
    m0col(m0) m1col(m1) k(2) s(4) method(exact)
matrix D_overlap = e(diagnostics)
// Target2 has two distinct equality blocks: propensities .20/.40 and
// .30/.35. Each block shares a 1/28 difference between adjacent raw rank
// weights. Arm masses are unchanged because each block lies within one arm;
// averaging raises mean propensity by ((.40-.20)+(.35-.30))/56=1/224.
matrix O_r = ( ///
    1,.5,.5,.5,1,.450357142857143,.35,.6 \ ///
    1,.285714285714286,.714285714285714,.285714285714286,1,(.4575+1/224),.2,.7 \ ///
    1,.535714285714286,.464285714285714,.535714285714286,1,.615,.4,.8)
mata: st_numscalar("_diag_r_overlap", max(abs(st_matrix("D_overlap")[,12..19] - st_matrix("O_r"))))
test_assert "[R parity] overlap diagnostics match RdistNN" scalar(_diag_r_overlap) < 1e-10

test_section "cate_dnn: [diagnostics] batch and grid default stability"

quietly cate_dnn y, xcols(x1 x2) x0matrix(X0_diag) tcol(D) ecol(e) ///
    m0col(m0) m1col(m1) k(1) tuningmode(grid)
matrix D_grid_default = e(diagnostics)
local s_default = e(s)
quietly cate_dnn y, xcols(x1 x2) x0matrix(X0_diag) tcol(D) ecol(e) ///
    m0col(m0) m1col(m1) k(1) s(`s_default')
matrix D_scalar_default = e(diagnostics)
mata: st_numscalar("_diag_grid_default", max(abs(st_matrix("D_grid_default") - st_matrix("D_scalar_default"))))
test_assert "[diagnostics] grid target diagnostics stay at default candidate" scalar(_diag_grid_default) < 1e-12

test_section "cate_dnn: [diagnostics] tied physical radii"

clear
input double x1 x2 rho_tie
-1  0  1
 1  0  2
 0 -1  3
 0  1  4
end
quietly cate_dnn rho_tie, xcols(x1 x2) x0(0 0) ///
    k(1) s(1) method(exact) pseudooutcome
matrix D_tie_radius = e(diagnostics)
test_assert "[diagnostics] tied 90-percent physical radius uses common Euclidean distance" ///
    abs(D_tie_radius[1,20] - 1) < 1e-12
test_assert "[diagnostics] score contribution diagnostics remain available for pseudo-outcomes" ///
    D_tie_radius[1,24] == 1 & D_tie_radius[1,25] < .

test_section "cate_dnn: [diagnostics] near-zero and zero arm mass"

clear
set obs 10
gen double x1 = _n
gen double y_arm = _n / 10
gen byte D_arm = (_n == 6)
gen double e_arm = .5
gen double m0_arm = 0
gen double m1_arm = 1
quietly cate_dnn y_arm, xcols(x1) x0(0) tcol(D_arm) ecol(e_arm) ///
    m0col(m0_arm) m1col(m1_arm) k(1) s(5)
matrix D_arm_edge = e(diagnostics)
test_assert "[diagnostics] tiny positive treated mass has normalized arm ESS one" ///
    D_arm_edge[1,13] > 0 & D_arm_edge[1,13] < .01 & ///
    abs(D_arm_edge[1,22] - 1) < 1e-12

replace D_arm = (_n == 10)
quietly cate_dnn y_arm, xcols(x1) x0(0) tcol(D_arm) ecol(e_arm) ///
    m0col(m0_arm) m1col(m1_arm) k(1) s(5)
matrix D_arm_zero = e(diagnostics)
test_assert "[diagnostics] zero treated local mass posts unavailable treated ESS" ///
    D_arm_zero[1,21] == 1 & D_arm_zero[1,13] == 0 & missing(D_arm_zero[1,22])

test_section "cate_dnn: [diagnostics] concentrated score contributions"

clear
set obs 20
gen double x1 = _n
gen double rho_concentrated = cond(_n == 1, 1e6, 0)
quietly cate_dnn rho_concentrated, xcols(x1) x0(0) ///
    k(1) s(5) method(exact) pseudooutcome
matrix D_score_concentrated = e(diagnostics)
test_assert "[diagnostics] top-five absolute score share detects concentration" ///
    D_score_concentrated[1,25] > .85
test_assert "[diagnostics] score share remains distinct from top-five raw-weight share" ///
    abs(D_score_concentrated[1,25] - D_score_concentrated[1,4]) > .01


test_section "cate_dnn: [diagnostics] underflowed within-arm squares"
// Literal powers of two give independent ESS (a+b)^2/(a^2+b^2)=9/5.
// The positive subnormal square fixture also loses the smaller square entirely.
mata {
    _essA=(0\1\1)
    _essSignals=J(3,2,0)
    _essDiag=J(1,25,.)
    _essDiag[1,1]=1
    for(_essPow=537;_essPow<=600;_essPow=_essPow+63) {
        _essW=(1\rdistnn_power2(-_essPow)\rdistnn_power2(-_essPow-1))
        for(_essFlip=0;_essFlip<=1;_essFlip++) {
            _essOut=rdistnn_clate_diagnostics(_essDiag,_essW,_essSignals,0,
                abs(_essA:-_essFlip),J(0,0,.))
            st_numscalar("_essLiteral"+strofreal(_essPow)+strofreal(_essFlip),
                _essOut[1,22+_essFlip])
            _essOut=rdistnn_clate_diagnostics(_essDiag,(1\.5\.25),_essSignals,0,
                abs(_essA:-_essFlip),J(0,0,.))
            st_numscalar("_essNormal"+strofreal(_essPow)+strofreal(_essFlip),
                _essOut[1,22+_essFlip])
            _essOut=rdistnn_clate_diagnostics(_essDiag,(1\_essW[2]\0),_essSignals,0,
                abs(_essA:-_essFlip),J(0,0,.))
            st_numscalar("_essSingleton"+strofreal(_essPow)+strofreal(_essFlip),
                _essOut[1,22+_essFlip])
        }
    }
    _essOut=rdistnn_clate_diagnostics(_essDiag,(1\0\0),_essSignals,0,
        _essA,J(0,0,.))
    st_numscalar("_essLiteralZero",_essOut[1,22])
    _essW=(.5\.25\.125)
    _essOut=rdistnn_clate_diagnostics(_essDiag,_essW,_essSignals,0,
        _essA,J(0,0,.))
    st_numscalar("_essOrdinarySame",_essOut[1,22]==quadsum(_essW[2..3])^2/quadsum(_essW[2..3]:^2))
}
foreach power in 537 600 {
    foreach flip in 0 1 {
        test_assert "[diagnostics] binary arm `flip' square exponent `power' recovers independent ESS 9/5" ///
            abs(scalar(_essLiteral`power'`flip')-9/5)<1e-14
        test_assert "[diagnostics] arm `flip' exponent `power' common scaling preserves ESS" ///
            abs(scalar(_essLiteral`power'`flip')-scalar(_essNormal`power'`flip'))<1e-14
        test_assert "[diagnostics] arm `flip' exponent `power' positive singleton has ESS one" ///
            scalar(_essSingleton`power'`flip')==1
    }
}
test_assert "[diagnostics] true zero-mass arm retains missing ESS" missing(scalar(_essLiteralZero))
test_assert "[diagnostics] ordinary arm ESS preserves established expression exactly" scalar(_essOrdinarySame)
