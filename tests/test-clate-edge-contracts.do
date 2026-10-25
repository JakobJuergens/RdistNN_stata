// Seeded R numerical audit, with full paired and cross-target covariance.
// Generate fixtures with Rscript tools/generate_r_clate_audit_fixtures.R.
// The generator itself validates independent subset/density-integral oracles.
version 17
do "tests/test_framework.do"
include "tests/fixtures/r_clate_audit_fixtures.do"

capture mata: mata drop ca_matrix_equal()
mata:
real scalar ca_matrix_equal(real matrix actual, real matrix expected,
    real scalar tol)
{
    real scalar scale
    real colvector a,b,ok
    if (rows(actual)!=rows(expected) | cols(actual)!=cols(expected)) return(0)
    a=vec(actual); b=vec(expected)
    if (any((a:>=.):!=(b:>=.))) return(0)
    ok=(b:<.)
    if (!any(ok)) return(1)
    a=select(a,ok); b=select(b,ok)
    scale=max(abs(b))
    if (scale==0) return(max(abs(a))<=1e-25)
    return(max(abs(a-b))/scale<=tol)
}
end

test_section "CLATE edge contracts: R goldens and independent paired-delete oracles"
forvalues case=1/`ca_cases' {
    clear
    svmat double CA_input`case'
    local p `ca_p`case''
    local xvars ""
    forvalues j=1/`p' {
        rename CA_input`case'`j' x`j'
        local xvars `xvars' x`j'
    }
    local next = `p'+1
    if `ca_nuisance`case'' {
        foreach name in Y D A e my0 my1 md0 md1 {
            rename CA_input`case'`next' `name'
            local ++next
        }
        local route "dcol(D) acol(A) ecol(e) my0col(my0) my1col(my1) md0col(md0) md1col(md1)"
        local outcome Y
    }
    else {
        rename CA_input`case'`next' py
        local ++next
        rename CA_input`case'`next' pd
        local route "pseudooutcome denominator(pd)"
        local outcome py
    }
    capture noisily clate_dnn `outcome', xcols(`xvars') ///
        x0matrix(CA_targets`case') `route' `ca_options`case''
    local fitrc = _rc
    test_assert "R audit case `case' succeeds: `ca_label`case''" `fitrc'==0
    if `fitrc' continue
    foreach pair in "estimates theta" "components components" ///
        "theta_vcov V" "component_vcov CV" "se_targets se" "status status" {
        local actual : word 1 of `pair'
        local expected : word 2 of `pair'
        mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(`actual')"),st_matrix("CA_`expected'`case'"),2e-8))
        test_assert "R audit case `case': full `actual' values and missing mask" scalar(ca_equal)==1
    }
}

test_section "shared singular Mahalanobis: independent geometry and reductions"
clear
svmat double CA_input12
rename (CA_input121 CA_input122 CA_input123 CA_input124) (x1 x2 py pd)
// X=(u,2u). Its Moore-Penrose covariance inverse is the literal rank-one
// formula below; an off-span target (0,1) projects to u=.4. No package
// inverse/rank/weight helper is used to construct these expected values.
weights_dnn x1 x2, x0(0 1) s(3) distance(mahalanobis) generate(ca_singular_w)
mata {
    ca_u=st_data(.,"x1")
    ca_uvar=quadsum((ca_u:-mean(ca_u)):^2)/7
    ca_mp=(1,2\2,4)/(25*ca_uvar)
    ca_expected_distance=abs(ca_u:-.4)/sqrt(ca_uvar)
    ca_expected_rank=(7\5\3\1\2\4\6\8)
    ca_expected_weight=(0\3\10\21\15\6\1\0)/56
    ca_geometry=rdistnn_mahal_geometry(st_data(.,"x1 x2"))
    ca_factor=ca_geometry[|3,1 \ rows(ca_geometry),cols(ca_geometry)|]/ca_geometry[1,1]
    st_numscalar("ca_equal",ca_matrix_equal(ca_factor*ca_factor',ca_mp,2e-12))
    st_numscalar("ca_distance_ok",ca_matrix_equal(rdistnn_distance_matrix(st_data(.,"x1 x2"),(0,1),"mahalanobis"),ca_expected_distance,2e-12))
    st_numscalar("ca_rank_ok",all(rdistnn_localization_ranks(rdistnn_distance_matrix(st_data(.,"x1 x2"),(0,1),"mahalanobis")):==ca_expected_rank))
    st_numscalar("ca_weight_ok",ca_matrix_equal(st_data(.,"ca_singular_w"),ca_expected_weight,2e-12))
}
test_assert "singular whitening factor Gram matches Moore-Penrose rank-one formula" scalar(ca_equal)==1
test_assert "off-span Mahalanobis distances match projected-coordinate formula" scalar(ca_distance_ok)==1
test_assert "off-span target rank order matches literal oracle" scalar(ca_rank_ok)==1
test_assert "standalone Mahalanobis weight command matches literal subset weights" scalar(ca_weight_ok)==1
np_dnn py, xcols(x1 x2) x0matrix(CA_targets12) s(3) distance(mahalanobis)
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(estimates)"),st_matrix("CA_theta12"),2e-12) & ca_matrix_equal(st_matrix("e(V)"),st_matrix("CA_V12"),2e-11))
test_assert "off-span singular NP estimates and full covariance match R ratio with unit denominator" scalar(ca_equal)==1
cate_dnn py, xcols(x1 x2) x0matrix(CA_targets12) s(3) distance(mahalanobis) pseudooutcome
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(estimates)"),st_matrix("CA_theta12"),2e-12) & ca_matrix_equal(st_matrix("e(V)"),st_matrix("CA_V12"),2e-11))
test_assert "off-span singular CATE estimates and full covariance match R ratio with unit denominator" scalar(ca_equal)==1

test_section "CLATE edge contracts: affine signal identities and target permutation"
clear
svmat double CA_input1
rename (CA_input11 CA_input12 CA_input13 CA_input14) (x1 x2 py pd)
clate_dnn py, denominator(pd) pseudooutcome xcols(x1 x2) ///
    x0matrix(CA_targets1) k(2) s(4)
matrix CA_baseB = e(estimates)
matrix CA_baseV = e(theta_vcov)
matrix CA_baseCV = e(component_vcov)
matrix CA_baseC = e(components)
gen double pa = -2*py+3*pd
clate_dnn pa, denominator(pd) pseudooutcome xcols(x1 x2) ///
    x0matrix(CA_targets1) k(2) s(4)
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(estimates)"),-2*st_matrix("CA_baseB"):+3,2e-12))
test_assert "affine numerator shifts target ratios" scalar(ca_equal)==1
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(theta_vcov)"),4*st_matrix("CA_baseV"),2e-11))
test_assert "affine numerator scales full ratio covariance" scalar(ca_equal)==1
mata {
    ca_T=I(8)
    for(ca_j=1;ca_j<=4;ca_j++) ca_T[(2*ca_j-1)..(2*ca_j),(2*ca_j-1)..(2*ca_j)]=(-2,3\0,1)
    st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(component_vcov)"),ca_T*st_matrix("CA_baseCV")*ca_T',2e-11))
}
test_assert "affine transformation preserves paired cross-component covariance" scalar(ca_equal)==1
gen double ps = -1e5*py
gen double ds = -1e5*pd
clate_dnn ps, denominator(ds) pseudooutcome xcols(x1 x2) ///
    x0matrix(CA_targets1) k(2) s(4)
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(estimates)"),st_matrix("CA_baseB"),2e-12) & ca_matrix_equal(st_matrix("e(theta_vcov)"),st_matrix("CA_baseV"),2e-11))
test_assert "common negative scaling preserves full ratio and covariance" scalar(ca_equal)==1
mata: st_numscalar("ca_equal",all(st_matrix("e(status)")[,2]:==1))
test_assert "common negative scaling updates denominator signs" scalar(ca_equal)==1
mata: st_matrix("CA_permuted",st_matrix("CA_targets1")[(4,2,1,3),])
clate_dnn py, denominator(pd) pseudooutcome xcols(x1 x2) ///
    x0matrix(CA_permuted) k(2) s(4)
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(estimates)"),st_matrix("CA_baseB")[,(4,2,1,3)],2e-12) & ca_matrix_equal(st_matrix("e(theta_vcov)"),st_matrix("CA_baseV")[(4,2,1,3),(4,2,1,3)],2e-11))
test_assert "target permutation carries full ratio covariance" scalar(ca_equal)==1
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(component_vcov)"),st_matrix("CA_baseCV")[(7,8,3,4,1,2,5,6),(7,8,3,4,1,2,5,6)],2e-11))
test_assert "target permutation carries interleaved full component covariance" scalar(ca_equal)==1
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(theta_vcov)")[,1],st_matrix("e(theta_vcov)")[,3],2e-12))
test_assert "duplicate targets retain identical covariance columns" scalar(ca_equal)==1

test_section "CLATE edge contracts: grid feasibility default and complete status masks"
clate_dnn py, denominator(pd) pseudooutcome xcols(x1 x2) ///
    x0matrix(CA_targets1) k(2) tuningmode(grid) ///
    sgrid(1 4 4 7 8 9) sdefault(4)
matrix CA_gridB = e(tuning_estimates)
matrix CA_gridSE = e(tuning_se)
matrix CA_gridY = e(tuning_reduced_form)
matrix CA_gridD = e(tuning_first_stage)
matrix CA_gridStatus = e(tuning_target_status)
matrix CA_gridMeta = e(tuning_candidates)
mata: st_numscalar("ca_equal",all(st_matrix("CA_gridB")[(1,5,6),]:>=.) & all(st_matrix("CA_gridSE")[(1,5,6),]:>=.) & all(st_matrix("CA_gridY")[(1,5,6),]:>=.) & all(st_matrix("CA_gridD")[(1,5,6),]:>=.) & all(st_matrix("CA_gridStatus")[(1..4,17..24),]:>=.))
test_assert "infeasible grid candidates retain every unavailable target and status" scalar(ca_equal)==1
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("CA_gridB")[2,],st_matrix("CA_theta1"),2e-10) & ca_matrix_equal(st_matrix("CA_gridSE")[2,],st_matrix("CA_se1"),2e-9) & ca_matrix_equal(st_matrix("CA_gridB")[2,],st_matrix("CA_gridB")[3,],2e-12))
test_assert "default and duplicate grid candidates match independent R fixture" scalar(ca_equal)==1
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(theta_vcov)"),st_matrix("CA_V1"),2e-9) & ca_matrix_equal(st_matrix("e(component_vcov)"),st_matrix("CA_CV1"),2e-9))
test_assert "grid active fit preserves full independent paired covariance" scalar(ca_equal)==1
test_assert "grid default ID and duplicate cache ID retained" e(tuning_default_candidate)==2 & el(CA_gridMeta,2,2)==el(CA_gridMeta,3,2)
clate_dnn py, denominator(pd) pseudooutcome xcols(x1 x2) ///
    x0matrix(CA_targets1) k(2) tuningmode(grid) sgrid(2 4 8) sdefault(8) nose
matrix CA_noseMeta = e(tuning_candidates)
matrix CA_noseStatus = e(tuning_target_status)
test_assert "nose grid permits s=n default" e(s)==8 & el(CA_noseMeta,3,4)==1
mata: st_numscalar("ca_equal",all(st_matrix("e(tuning_se)"):>=.) & all(st_matrix("CA_noseStatus")[,6]:==1) & all(st_matrix("CA_noseStatus")[,7]:==0))
test_assert "nose grid consistently marks covariance uncomputed" scalar(ca_equal)==1

test_section "CLATE edge contracts: nuisance relabeling and sample boundary validation"
clear
svmat double CA_input26
rename (CA_input261 CA_input262 CA_input263 CA_input264 CA_input265 ///
    CA_input266 CA_input267 CA_input268 CA_input269 CA_input2610) ///
    (x1 x2 Y D A e my0 my1 md0 md1)
local nuis "xcols(x1 x2) x0matrix(CA_targets26) dcol(D) acol(A) ecol(e) my0col(my0) my1col(my1) md0col(md0) md1col(md1) k(2) s(4) clip(.05) distance(mahalanobis)"
clate_dnn Y, `nuis'
matrix CA_ivB = e(estimates)
matrix CA_ivV = e(theta_vcov)
matrix CA_ivC = e(components)
gen byte Af = 1-A
gen double ef = 1-e
clate_dnn Y, xcols(x1 x2) x0matrix(CA_targets26) dcol(D) acol(Af) ecol(ef) ///
    my0col(my1) my1col(my0) md0col(md1) md1col(md0) ///
    k(2) s(4) clip(.05) distance(mahalanobis)
mata: st_numscalar("ca_equal",ca_matrix_equal(st_matrix("e(estimates)"),st_matrix("CA_ivB"),2e-11) & ca_matrix_equal(st_matrix("e(theta_vcov)"),st_matrix("CA_ivV"),2e-10) & ca_matrix_equal(st_matrix("e(components)"),-st_matrix("CA_ivC"),2e-11))
test_assert "relabeling instrument arms negates both components and preserves ratio covariance" scalar(ca_equal)==1
foreach var in e md0 md1 {
    gen double save_`var' = `var'
    replace `var' = 1.2 in 1
    capture clate_dnn Y if _n>1, `nuis'
    local restricted_rc = _rc
    test_assert "excluded invalid `var' prediction does not enter fit" `restricted_rc'==0
    capture clate_dnn Y, `nuis'
    local included_rc = _rc
    test_assert "included invalid `var' prediction fails before clipping" `included_rc'==198
    replace `var' = save_`var'
}
replace A = 2 in 1
capture clate_dnn Y, `nuis'
local binary_rc = _rc
test_assert "nonbinary instrument rejected" `binary_rc'==198
replace A = 0 in 1
replace Y = . in 1
clate_dnn Y, `nuis'
test_assert "missing outcome removed from estimation sample" e(N)==9
gen byte ca_esample = e(sample)
test_assert "outcome missing row absent from e(sample)" ca_esample[1]==0

// Public invalid targets must fail cleanly rather than manufacture zero ratios.
matrix CA_badtarget = (0,. \ 1,1)
capture clate_dnn Y, xcols(x1 x2) x0matrix(CA_badtarget) ///
    dcol(D) acol(A) ecol(e) my0col(my0) my1col(my1) md0col(md0) md1col(md1) s(3)
local target_rc = _rc
test_assert "missing target coordinate rejected" `target_rc'==198


test_section "CLATE edge contracts: relative denominator status at ordinary cutoffs"
clear
set obs 6
gen double cr_x=_n
gen double cr_num=1
gen double cr_den=1
replace cr_den=4 in 6
// Exact s=n weights select row 1. The independently known relative magnitude
// is cr_den[1]/4; neighboring represented values must straddle the same cutoff.
foreach sign in 1 -1 {
    foreach side in -1 0 1 {
        replace cr_den=`sign'*(1+`side'*2^-50) in 1
        quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
            xcols(cr_x) x0(0) s(6) nose denominatortol(.25)
        local cr_want=(`side'<=0)
        test_assert "ordinary signed cutoff neighbor `sign' `side' preserves exact decision" ///
            el(e(status),1,3)==`cr_want' & el(e(status),1,1)==0 & ///
            el(e(status),1,2)==(`sign'<0) & ///
            el(e(estimates),1,1)==1/cr_den[1]
    }
}
replace cr_den=0 in 1
quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
    xcols(cr_x) x0(0) s(6) nose denominatortol(.25)
test_assert "zero denominator remains undefined and is not a relative near-zero flag" ///
    el(e(status),1,1)==1 & el(e(status),1,3)==0 & missing(el(e(estimates),1,1))

// A positive binary common scale makes the independent reference safe without
// taking abs() of the original large negative value. Row 6 has zero weight.
test_section "CLATE edge contracts: finite negative denominator magnitudes"
replace cr_den=1
replace cr_den=-8e307-8e307 in 6
scalar cr_small_relative=(1/4)/(4e307)
test_assert "independent inactive-signal relative reference is positive and below tolerance" ///
    cr_small_relative>0 & cr_small_relative<1e-8
foreach common in 1 4 {
    replace cr_num=1/`common'
    replace cr_den=1/`common'
    replace cr_den=(-8e307-8e307)/`common' in 6
    quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
        xcols(cr_x) x0(0) s(6) nose denominatortol(1e-8)
    test_assert "inactive huge negative denominator defines all-input scale /`common'" ///
        el(e(estimates),1,1)==1 & el(e(components),1,2)==1/`common' & ///
        el(e(status),1,1)==0 & el(e(status),1,3)==1
    quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
        xcols(cr_x) x0(0) s(6) nose denominatortol(0)
    test_assert "zero tolerance disables exceptional relative flag /`common'" ///
        el(e(estimates),1,1)==1 & el(e(status),1,3)==0
}
// This time the fitted local denominator itself has unsafe positive magnitude.
// Its magnitude relative to the all-input maximum is exactly one.
replace cr_den=-8e307-8e307
replace cr_num=cr_den
foreach common in 1 4 {
    foreach cutoff in .5 1 {
        quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
            xcols(cr_x) x0(0) s(6) nose denominatortol(`cutoff')
        test_assert "huge negative local denominator compares relative one at `cutoff' /`common'" ///
            el(e(estimates),1,1)==1 & el(e(status),1,1)==0 & ///
            el(e(status),1,2)==1 & el(e(status),1,3)==(`cutoff'==1)
    }
    replace cr_num=cr_num/4
    replace cr_den=cr_den/4
}

// Only the descriptive flag may change when tolerance changes. The existing
// independent four-row oracle also protects finite paired component covariance.
clear
set obs 4
gen double cr_x=_n
gen double cr_num=cond(_n==1,2,cond(_n==2,1,cond(_n==3,4,3)))
gen double cr_den=cond(_n==1,1,cond(_n==2,-1,cond(_n==3,2,1)))
quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
    xcols(cr_x) x0(0) s(2) denominatortol(0)
matrix CR_B=e(estimates)
matrix CR_C=e(components)
matrix CR_V=e(theta_vcov)
matrix CR_CV=e(component_vcov)
quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
    xcols(cr_x) x0(0) s(2) denominatortol(1)
mata: st_numscalar("cr_same",all(st_matrix("CR_B"):==st_matrix("e(estimates)")) & all(st_matrix("CR_C"):==st_matrix("e(components)")) & all(st_matrix("CR_V"):==st_matrix("e(theta_vcov)")) & all(st_matrix("CR_CV"):==st_matrix("e(component_vcov)")))
test_assert "relative flag preserves estimates and all paired component covariance bytes" ///
    cr_same & el(e(status),1,3)==1 & el(CR_B,1,1)==4 & ///
    abs(el(CR_V,1,1)-34)<1e-10 & abs(el(CR_CV,1,2)-.5)<1e-12

// Each route evaluates the known relative magnitudes (tiny, one), with original
// signals retained and duplicate/reordered target identities kept explicitly.
test_section "CLATE edge contracts: relative status survives prediction candidates and cold restore"
clear
set obs 6
gen double cr_x=_n
gen double cr_den=1
replace cr_den=-8e307-8e307 in 6
gen double cr_num=cr_den
matrix CR_TARGETS=(0\7\0)
matrix CR_REORDER=(7\0\7\0)
foreach common in 1 4 {
    quietly clate_dnn cr_num, denominator(cr_den) pseudooutcome ///
        xcols(cr_x) x0matrix(CR_TARGETS) tuningmode(grid) ///
        sgrid(5 6 6) sdefault(6) nose denominatortol(1e-8)
    mata: st_numscalar("cr_status_ok",all(st_matrix("e(status)")[,3]:==(1\0\1)) & all(st_matrix("e(estimates)"):==1) & all(st_matrix("e(tuning_target_status)")[,5]:==(1\0\1\1\0\1\1\0\1)))
    test_assert "original and duplicate grid candidates use independent relative references /`common'" cr_status_ok
    tempname cr_parent
    estimates store `cr_parent'
    quietly rdistnn_candidate, candidate(3) post
    mata: st_numscalar("cr_status_ok",all(st_matrix("e(status)")[,3]:==(1\0\1)) & all(st_matrix("e(estimates)"):==1))
    test_assert "duplicate candidate retains corrected denominator statuses /`common'" cr_status_ok
    quietly estimates restore `cr_parent'
    quietly rdistnn_predict, x0matrix(CR_REORDER) post
    mata: st_numscalar("cr_status_ok",all(st_matrix("e(status)")[,3]:==(0\1\0\1)) & all(st_matrix("e(estimates)"):==1) & all(st_matrix("e(status)")[,2]:==(1\0\1\0)))
    test_assert "prediction reorders and repeats relative statuses with target identities /`common'" cr_status_ok
    quietly estimates restore `cr_parent'
    tempfile cr_saved
    quietly estimates save "`cr_saved'", replace
    preserve
    clear
    mata: mata clear
    quietly estimates use "`cr_saved'"
    quietly rdistnn_predict, x0matrix(CR_REORDER) post
    mata: st_numscalar("cr_status_ok",all(st_matrix("e(status)")[,3]:==(0\1\0\1)) & all(st_matrix("e(estimates)"):==1))
    test_assert "cold restored prediction retains all original denominator signals /`common'" cr_status_ok & _N==0
    restore
    estimates drop `cr_parent'
    replace cr_num=cr_num/4
    replace cr_den=cr_den/4
}
