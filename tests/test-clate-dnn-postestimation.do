// CLATE adapters: canonical target tables, paired influence, and shared plots.
version 17
capture mkdir "build"
capture mkdir "build/clate"
do "tests/test_framework.do"
test_section "CLATE postestimation: cold loading and full workflow"
clear
set obs 90
set seed 3219
gen double x = rnormal()
gen byte a = runiform() < .5
gen byte d = a*(runiform()<.8)
gen double ep = .5
gen double my0 = x
gen double my1 = x+1.2
gen double md0 = 0
gen double md1 = .8
gen double y = x+1.5*d+rnormal()
gen double sy = my1-my0+a*(y-my1)/ep-(1-a)*(y-my0)/(1-ep)
gen double sd = md1+a*(d-md1)/ep
matrix CLX = (-1 \ -.5 \ 0 \ .5 \ 1)
quietly clate_dnn y, x(x) x0matrix(CLX) dcol(d) acol(a) ecol(ep) ///
    my0col(my0) my1col(my1) md1col(md1) onesided k(2) tuningmode(grid)
estimates store clpost
matrix CLB = e(estimates)
capture noisily clate_dnn_plot, name(clcold) nodraw replace
local coldrc = _rc
capture matrix CLPOS = r(plot_data)
local colddatarc = _rc
test_assert "[clate cold] public plot loads shared builders" `coldrc' == 0
if `coldrc'==0 & `colddatarc'==0 {
    matrix CLSE = e(se_targets)
    mata: _clpos = st_matrix("CLPOS"); _clse = st_matrix("CLSE")'; _clb = st_matrix("CLB")'
    mata: st_numscalar("clpositive_error", max(abs((_clpos[,7]-_clse, _clpos[,8]-(_clb-invnormal(.975)*_clse), _clpos[,9]-(_clb+invnormal(.975)*_clse)))))
    mata: st_numscalar("clpositive_width", min(_clpos[,9]-_clpos[,8]))
    test_assert "[clate Wald] positive available SEs retain canonical confidence bounds" scalar(clpositive_error)<1e-12 & scalar(clpositive_width)>0
}
else test_assert "[clate Wald] positive available SE plot data returned" 0
capture noisily clate_dnn_diagnostics, targets(2 4)
local diagrc = _rc
test_assert "[clate diagnostics] target selection supported" `diagrc' == 0
if `diagrc'==0 {
    matrix CLDIAG = r(diagnostics)
    matrix CLCOMP = r(components)
    matrix CLSTAT = r(status)
    test_assert "[clate diagnostics] selected target identities retained" CLDIAG[1,1]==2 & CLDIAG[2,1]==4
    test_assert "[clate diagnostics] selected component table includes first-stage uncertainty" rowsof(CLCOMP)==2 & colsof(CLCOMP)==9
    test_assert "[clate diagnostics] selected status table retains five codes" rowsof(CLSTAT)==2 & colsof(CLSTAT)==5
    matrix CLV = e(component_vcov)
    test_assert "[clate diagnostics] denominator SE follows paired component covariance" abs(CLCOMP[1,5]-sqrt(CLV[4,4]))<1e-12
    test_assert "[clate diagnostics] signed first-stage statistic follows component SE" abs(CLCOMP[1,8]-CLCOMP[1,3]/CLCOMP[1,5])<1e-12
    local dnames : colnames CLDIAG
    test_assert "[clate diagnostics] instrument arm names exposed" strpos("`dnames'", "instrument_1_mass")>0
}
capture noisily clate_dnn_diagnostics, candidates(1 3)
local diagrc = _rc
test_assert "[clate diagnostics] candidate selection supported" `diagrc'==0
if `diagrc'==0 {
    local crows = r(n_rows)
    test_assert "[clate diagnostics] two candidate rows" `crows'==2
}
capture noisily clate_dnn_diagnostics, allcandidates
test_assert "[clate diagnostics] all candidate diagnostics supported" _rc==0
tempfile cldiag
capture noisily clate_dnn_diagnostics, targets(1 5) saving("`cldiag'") replace
local saverc = _rc
test_assert "[clate diagnostics] saving selected targets works" `saverc'==0
if `saverc'==0 {
    preserve
    use "`cldiag'", clear
    test_assert "[clate diagnostics] saved table has component and status columns" _N==2 & !missing(first_stage[1]) & ratio_status[1]==0
    restore
}
foreach kind in estimates diagnostics tuning overlap weights influence dashboard {
    capture noisily clate_dnn_plot, type(`kind') name(cl_`kind') replace nodraw
    test_assert "[clate plot] `kind' builds" _rc==0
}
foreach style in auto profile pointrange forest {
    capture noisily clate_dnn_plot, style(`style') targets(1 3 5) name(clstyle) replace nodraw
    local stylerc = _rc
    if `stylerc'==0 & "`style'"=="forest" matrix CLSELECT = r(plot_data)
    test_assert "[clate plot] estimate style `style' builds" `stylerc'==0
}
test_assert "[clate plot data] estimates return selected target rows and identities" rowsof(CLSELECT)==3 & CLSELECT[1,5]==1 & CLSELECT[2,5]==3 & CLSELECT[3,5]==5
capture noisily clate_dnn_plot, type(tuning) targets(2 4) candidates(1 3) name(cltuneselect) replace nodraw
local tuneselectrc = _rc
capture matrix CLTSELECT = r(plot_data)
local tunedatarc = _rc
test_assert "[clate plot data] tuning returns selected target and candidate identities" `tuneselectrc'==0 & `tunedatarc'==0 & rowsof(CLTSELECT)==4 & CLTSELECT[1,1]==1 & CLTSELECT[3,1]==3 & CLTSELECT[1,5]==2 & CLTSELECT[2,5]==4
foreach view in mass armess propensity {
    capture noisily clate_dnn_plot, type(overlap) view(`view') targets(1 3) name(cloverlap) replace nodraw
    test_assert "[clate plot] instrument overlap `view' builds" _rc==0
}
foreach view in profile cumulative kernel {
    local select ""
    if "`view'"=="kernel" local select "targets(1 3)"
    capture noisily clate_dnn_plot, type(weights) view(`view') `select' name(clweights) replace nodraw
    test_assert "[clate plot] weights `view' builds" _rc==0
}
foreach kind in estimates diagnostics tuning weights {
    capture noisily clate_dnn_plot, type(`kind') candidates(1 3) name(clcandidates) replace nodraw
    test_assert "[clate plot] `kind' candidate selection works" _rc==0
}
capture noisily clate_dnn_plot, type(influence) targets(2 4) name(clinfluence) nlabel(0) replace nodraw
local infrc = _rc
local infn = r(n_plotted)
if `infrc'==0 matrix CLIF = r(plot_data)
test_assert "[clate influence] selected targets retained" `infrc'==0 & `infn'==2
if `infrc'==0 {
    mata {
        _cli = st_matrix("CLIF")
        _clb = st_matrix("CLB")
        _clsy = st_data(.,"sy")
        _clsd = st_data(.,"sd")
        _cldiff = 0
        for (_clq=1;_clq<=rows(_cli);_clq++) {
            _clobs = _cli[_clq,2]
            _clt = _cli[_clq,1]
            _clref = _clsy[_clobs] - _clb[1,_clt]*_clsd[_clobs]
            _cldiff = max((_cldiff, abs(_cli[_clq,6]-_clref),
                abs(_cli[_clq,7]-_cli[_clq,5]*_clref),
                abs(_cli[_clq,8]-_clsy[_clobs]), abs(_cli[_clq,9]-_clsd[_clobs])))
        }
        st_numscalar("clpaired_error",_cldiff)
    }
    test_assert "[clate influence] reconstructed paired contributions match independent raw signals" scalar(clpaired_error)<1e-12
}
// A unit denominator cannot detect accidentally discarding the second signal.
// Verify reconstructed contribution data using target-specific centering.
quietly clate_dnn sy, denominator(sd) pseudooutcome x(x) x0matrix(CLX) s(12)
capture noisily clate_dnn_plot, type(influence) name(clpaired) replace nodraw
test_assert "[clate influence] precomputed paired signals work" _rc==0
quietly clate_dnn y, x(x) x0matrix(CLX) dcol(d) acol(a) ecol(ep) ///
    my0col(my0) my1col(my1) md1col(md1) onesided k(2) tuningmode(grid)
foreach ext in png svg gph {
    capture noisily clate_dnn_plot, name(clexport) replace nodraw saving("build/clate/post-estimates.`ext'")
    local exprc = _rc
    capture confirm file "build/clate/post-estimates.`ext'"
    test_assert "[clate export] `ext' artifact exists" `exprc'==0 & _rc==0
}
tempname clsvg
file open `clsvg' using "build/clate/post-estimates.svg", read text
file read `clsvg' clline
local clpolygon 0
while r(eof)==0 {
    if strpos(`"`clline'"',"<path id=")>0 & strpos(`"`clline'"'," L")>0 local clpolygon 1
    file read `clsvg' clline
}
file close `clsvg'
test_assert "[clate export] confidence band is a rendered polygon" `clpolygon'==1
capture noisily clate_dnn_plot, type(dashboard) name(cldash) replace nodraw saving("build/clate/post-dashboard.png")
test_assert "[clate export] dashboard artifact exists" _rc==0
quietly estimates restore clpost
matrix CLAFTER = e(estimates)
mata: st_numscalar("cldiff", max(abs(st_matrix("CLAFTER")-st_matrix("CLB"))))
test_assert "[clate restore] estimation matrices and original data retained" scalar(cldiff)<1e-12 & _N==90
capture noisily clate_dnn
 test_assert "[clate replay] stored result replays after plots" _rc==0
preserve
replace y = y + 100 in 1
capture noisily clate_dnn_plot, type(influence) name(clmutated) replace nodraw
test_assert "[clate signature] outcome mutation rejects influence" _rc==498
capture noisily clate_dnn_plot, type(weights) view(kernel) name(clmutated) replace nodraw
test_assert "[clate signature] outcome mutation rejects kernel" _rc==498
capture noisily clate_dnn_plot, type(estimates) name(clstored) replace nodraw
test_assert "[clate signature] stored estimate plots remain usable after data change" _rc==0
restore
capture noisily clate_dnn_plot, type(influence) targets(2) name(clrestored) replace nodraw
test_assert "[clate signature] restoring original data restores influence" _rc==0
preserve
gsort -x
capture noisily clate_dnn_plot, type(influence) name(clsorted) replace nodraw
test_assert "[clate signature] changed data order rejects influence" _rc==498
restore
preserve
drop md1
capture noisily clate_dnn_plot, type(influence)
local stale = _rc
restore
test_assert "[clate influence] missing original nuisance rejected" `stale'==498
// An all-undefined fit has no e(b); plotting must use full canonical tables.
gen double zero = 0
quietly clate_dnn sy, denominator(zero) pseudooutcome x(x) x0matrix(CLX) s(12)
capture noisily clate_dnn_plot, name(clundefined) replace nodraw
test_assert "[clate undefined] estimates work without e(b)" _rc==0
capture noisily clate_dnn_diagnostics, targets(2 4)
test_assert "[clate undefined] diagnostic targets retained" _rc==0
capture noisily clate_dnn_plot, type(influence) targets(1)
test_assert "[clate undefined] influence reports undefined ratio" _rc==498
capture noisily clate_dnn_plot, type(overlap)
test_assert "[clate overlap] paired route reports unavailable instrument data" _rc==498
// Canonical target identities survive partial e(b) omission.
clear
set obs 4
gen double x = _n
gen double sy = 1
gen double sd = cond(_n<=2,0,1)
matrix CLPART = (0 \ 5)
quietly clate_dnn sy, denominator(sd) pseudooutcome x(x) x0matrix(CLPART) s(4) nose
capture noisily clate_dnn_plot, style(forest) name(clpartial) replace nodraw
local partialrc = _rc
local partialn = r(n_plotted)
test_assert "[clate partial] full-target plot does not index reduced e(b)" `partialrc'==0 & `partialn'==2
capture noisily clate_dnn_plot, type(influence) targets(2) name(clpartialinf) replace nodraw
test_assert "[clate partial] defined target influence uses canonical identity" _rc==0
capture noisily clate_dnn_plot, type(influence) targets(1)
test_assert "[clate partial] undefined target influence rejected" _rc==498
// Zero paired score variance retains a ratio/SE but no Wald interval.
clear
set obs 20
gen double x = _n
gen double sy = 2
gen double sd = 1
quietly clate_dnn sy, denominator(sd) pseudooutcome x(x) x0(10) tuningmode(grid)
capture noisily clate_dnn_plot, name(clwald) replace nodraw
local waldrc = _rc
capture matrix CWALD = r(plot_data)
local walddata = _rc
test_assert "[clate Wald] estimate points retained but unavailable CIs missing" `waldrc'==0 & `walddata'==0 & CWALD[1,6]==2 & missing(CWALD[1,8]) & missing(CWALD[1,9])
capture noisily clate_dnn_plot, type(tuning) name(clwaldtune) replace nodraw
local waldrc = _rc
capture matrix CWALDT = r(plot_data)
local walddata = _rc
test_assert "[clate Wald] tuning points retained but unavailable CIs missing" `waldrc'==0 & `walddata'==0 & CWALDT[1,6]==2 & missing(CWALDT[1,8]) & missing(CWALDT[1,9])
ereturn clear
capture noisily clate_dnn_plot
test_assert "[clate errors] plot requires CLATE estimates" _rc==301
capture noisily clate_dnn_diagnostics
test_assert "[clate errors] diagnostics require CLATE estimates" _rc==301
