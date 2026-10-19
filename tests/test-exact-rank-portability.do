version 17
do "tests/test_framework.do"
quietly do "src/_rdistnn_mata.ado"
test_section "exact rank-table portability"
clear
import delimited "tests/fixtures/exact_rank_tables.csv", clear varnames(1) asdouble
mata:
_fn_cases=uniqrows(st_data(.,("n","k","s")))
_fn_ok=1
for(_fn_i=1;_fn_i<=rows(_fn_cases);_fn_i++) {
    _fn_case=_fn_cases[_fn_i,]
    _fn_ix=selectindex(rowsum(st_data(.,("n","k","s")):==_fn_case):==3)
    _fn_w=rdistnn_exact_table(_fn_case[1],_fn_case[2],_fn_case[3])
    if(!all(_fn_w:==st_data(_fn_ix,"weight"))) _fn_ok=0
}
st_numscalar("_fn_exact",_fn_ok)
end
test_assert "Mata exact tables match actual R compiled and fallback doubles" scalar(_fn_exact)
clear
set obs 20
gen double fx=_n
gen double fy=.
mata: st_store(.,"fy",(17\2\13\9\6\10\4\19\3\12\7\18\20\1\16\11\14\15\5\8))
foreach engine in mata plugin {
    if "`engine'"=="plugin" global rdistnn_use_plugin 1
    else global rdistnn_use_plugin 0
    quietly np_dnn fy, xcols(fx) x0(0) s(2) nose
    matrix FQ=e(target_diagnostics)
    test_assert "`engine': end-to-end exact Q90 uses common represented table" FQ[1,17]==19
}
macro drop rdistnn_use_plugin
mata:
_fn_n=2000; _fn_k=_fn_n-2
_fn_w=rdistnn_exact_table(_fn_n,_fn_k,_fn_n-1)
_fn_oracle=(J(_fn_k,1,(_fn_n-1)/_fn_n/_fn_k)\1/_fn_n\0)
st_numscalar("_fn_sparse",max(abs(_fn_w-_fn_oracle))<3e-18)
end
test_assert "near-full subsets retain analytic weights with bounded support work" scalar(_fn_sparse)
