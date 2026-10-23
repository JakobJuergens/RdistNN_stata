// Run in a fresh Stata process per installation. This benchmark is evidence,
// not a CI timing threshold. Repetitions default to 3; use the same value for
// before/after. Full numerical matrices are retained outside package source.
// do tools/bench/bench_reporting.do "installation-root" "receipt-root" "before" 3
version 17
args installation output phase repetitions
if `"`installation'"'=="" | `"`output'"'=="" | `"`phase'"'=="" exit 198
if "`repetitions'"=="" local repetitions 3
set more off
set matsize 11000
local installation=subinstr("`installation'","\","/",.)
local output=subinstr("`output'","\","/",.)
adopath ++ "`installation'/src"
capture mkdir "`output'"
capture mkdir "`output'/`phase'"
quietly do "`installation'/src/_rdistnn_mata.ado"
// Instrument source copies in the receipt directory, leaving both installed
// trees untouched. The numerical engines still execute their real code.
capture mkdir "`installation'/build"
local shadow "`installation'/build/reporting-bench-instrument"
capture mkdir "`shadow'"
mata:
void bench_copy(string scalar source, string scalar destination, real scalar blocks)
{
    real scalar input,output,in_notice,depth
    string scalar line
    input=fopen(source,"r"); unlink(destination); output=fopen(destination,"w")
    in_notice=0; depth=0
    while((line=fget(input))!=J(0,0,"")) {
        if(blocks) {
            if(strpos(line,"void rdistnn_distance_warning(")) in_notice=1
            if(in_notice & strpos(line,"rdistnn_duplicate_rows(DX[,j])")) {
                line=subinstr(line,"rdistnn_duplicate_rows(DX[,j])","bench_notice_duplicates(DX[,j])")
            }
            line=subinstr(line,"n=rows(ranks); m=cols(ranks)","st_numscalar("+char(34)+"BENCH_GROUP_CALLS"+char(34)+",st_numscalar("+char(34)+"BENCH_GROUP_CALLS"+char(34)+")+1); n=rows(ranks); m=cols(ranks)")
            if(in_notice & line=="}") in_notice=0
        }
        else line=subinstr(line,"st_global("+char(34)+"rdistnn_last_distance_engine"+char(34)+","+char(34)+char(34)+")","st_numscalar("+char(34)+"BENCH_DISTANCE_CALLS"+char(34)+",st_numscalar("+char(34)+"BENCH_DISTANCE_CALLS"+char(34)+")+1); st_global("+char(34)+"rdistnn_last_distance_engine"+char(34)+","+char(34)+char(34)+")")
        fput(output,line)
    }
    fclose(input); fclose(output)
}
real scalar bench_notice_duplicates(real matrix X)
{
    st_numscalar("BENCH_NOTICE_SORTS",st_numscalar("BENCH_NOTICE_SORTS")+1)
    return(rdistnn_duplicate_rows(X))
}
bench_copy(st_local("installation")+"/src/_rdistnn_blocks.ado",st_local("shadow")+"/_rdistnn_blocks.ado",1)
bench_copy(st_local("installation")+"/src/_rdistnn_distance_mata.ado",st_local("shadow")+"/_rdistnn_distance_mata.ado",0)
end
adopath ++ "`shadow'"
tempname timing
postfile `timing' str8 phase long n int m byte tied str6 engine str8 output_mode str10 route byte repetition double seconds long group_calls distance_calls notice_sorts str12 actual_distance str12 actual_deletion using "`output'/`phase'/timings.dta", replace every(1)
foreach n in 64 2000 32000 {
    local m=cond(`n'==64,1,cond(`n'==2000,8,16))
    foreach tied in 0 1 {
        clear
        quietly set obs `n'
        generate double x1=cond(`tied',floor((_n-1)/4)+1,_n)
        generate double y=sin(_n/10)
        mata: BENCH_SEQ=(0::(`m'-1))/max((1,`m'-1)); BENCH_TARGETS=(.125:+.25*BENCH_SEQ); st_matrix("BENCH_TARGETS",BENCH_TARGETS)
        foreach engine in mata plugin {
            global rdistnn_distance_engine `engine'
            global rdistnn_jackknife_engine `engine'
            foreach output_mode in ordinary quiet {
                local prefix ""
                if "`output_mode'"=="quiet" local prefix "quietly"
                local routes "fit"
                if `n'==2000 local routes "fit grid dashboard"
                foreach route of local routes {
                    forvalues repetition=1/`repetitions' {
                        // Dashboard reuses the preceding grid's retained fit;
                        // it does not refit geometry outside the timed operation.
                        scalar BENCH_GROUP_CALLS=0
                        scalar BENCH_DISTANCE_CALLS=0
                        scalar BENCH_NOTICE_SORTS=0
                        timer clear 91
                        timer on 91
                        if "`route'"=="fit" `prefix' np_dnn y, xcols(x1) x0matrix(BENCH_TARGETS) s(8) nose
                        if "`route'"=="grid" `prefix' np_dnn y, xcols(x1) x0matrix(BENCH_TARGETS) tuningmode(grid) sgrid(2 3 4 5 6 7 8 9 10 11 12) sdefault(8) nose
                        if "`route'"=="dashboard" {
                            if `n'*`m'<=11000 `prefix' np_dnn_plot, type(dashboard) name(bench_reporting) replace nodraw
                            else {
                                // Stata/SE's graph mkmat transport caps rows at
                                // 11,000. Exercise the real three Mata panel
                                // preparations at larger geometries instead.
                                quietly findfile "_rdistnn_model_mata.ado"
                                quietly run "`r(fn)'"
                                capture program list _rdistnn_report_notices
                                if !_rc `prefix' _rdistnn_report_notices
                                `prefix' mata: BENCH_PANEL_W=rdistnn_model_weights((1..`m')); BENCH_PANEL_D=rdistnn_plot_target_diag(st_matrix("e(diagnostics)"),st_matrix("e(targets)")); BENCH_PANEL_E=rdistnn_plot_est_default(st_matrix("e(targets)"),st_matrix("e(estimates)")',st_matrix("e(se_targets)")',8,invnormal(.975))
                            }
                        }
                        timer off 91
                        local actual_distance "$rdistnn_last_distance_engine"
                        local actual_deletion "$rdistnn_last_jackknife_engine"
                        if inlist("`route'","fit","grid") & "`actual_distance'"!="`engine'" exit 498
                        local actual_deletion "none"
                        if "`route'"=="dashboard" local actual_distance "cached"
                        local groups=BENCH_GROUP_CALLS
                        local distances=BENCH_DISTANCE_CALLS
                        local notice_sorts=BENCH_NOTICE_SORTS
                        // Numerical receipts are outside timing and retain every
                        // weight/rank/group plus estimates/covariance/candidate rows.
                        // Intentional nose matches the reporting-cost R benchmark;
                        // deletion-engine parity has separate runtime/native gates.
                        if `repetition'==1 {
                            mata: BENCH_W=_np_W; BENCH_R=st_matrix("e(model_ranks)"); BENCH_G=st_matrix("e(model_groups)"); BENCH_D=st_matrix("e(model_distances)"); BENCH_B=st_matrix("e(estimates)"); BENCH_V=st_matrix("e(estimate_vcov)"); BENCH_C=BENCH_CE=BENCH_CD=J(0,0,.)
                            if inlist("`route'","grid","dashboard") {
                                mata: BENCH_C=st_matrix("e(tuning_candidates)"); BENCH_CE=st_matrix("e(tuning_estimates)"); BENCH_CD=st_matrix("e(tuning_diagnostics)")
                            }
                            mata: mata matsave "`output'/`phase'/`n'-`tied'-`engine'-`output_mode'-`route'.mmat" BENCH_W BENCH_R BENCH_G BENCH_D BENCH_B BENCH_V BENCH_C BENCH_CE BENCH_CD, replace
                        }
                        quietly timer list 91
                        post `timing' ("`phase'") (`n') (`m') (`tied') ("`engine'") ("`output_mode'") ("`route'") (`repetition') (r(t91)) (`groups') (`distances') (`notice_sorts') ("`actual_distance'") ("`actual_deletion'")
                    }
                }
            }
        }
    }
}
postclose `timing'
use "`output'/`phase'/timings.dta", clear
export delimited using "`output'/`phase'/timings.csv", replace
// The Windows launcher separately records process peak working-set bytes.
// Stata does not expose allocation traffic; no allocation count is invented.
display "REPORTING BENCHMARK COMPLETE: `phase'"
