do "tests/test_framework.do"
clear
set obs 30
gen double x = (_n-15)/4
gen double py = 1+sin(x)
gen double pd = .4+cos(x)/10
matrix CT = (-1 \ 0 \ 1)
test_section "clate_dnn: shared tuning grid and methods"
foreach method in exact asymp {
    foreach distance in l2 mahalanobis {
        clate_dnn py, denominator(pd) pseudooutcome xcols(x) ///
            x0matrix(CT) method(`method') distance(`distance') ///
            k(2) tuningmode(grid) sgrid(1 4 4 8 30 31) sdefault(4)
        matrix GB = e(tuning_estimates)
        matrix GS = e(tuning_se)
        matrix GY = e(tuning_reduced_form)
        matrix GD = e(tuning_first_stage)
        matrix GM = e(tuning_candidates)
        matrix B = e(estimates)
        test_assert "grid rows preserved `method'/`distance'" rowsof(GB)==6
        test_assert "infeasible low k retained `method'/`distance'" el(GM,1,4)==0
        test_assert "infeasible full order retained `method'/`distance'" el(GM,5,4)==0
        test_assert "duplicate order values equal `method'/`distance'" abs(el(GB,2,1)-el(GB,3,1))<1e-12
        test_assert "duplicate computed IDs equal `method'/`distance'" el(GM,2,2)==el(GM,3,2)
        test_assert "default candidate unchanged `method'/`distance'" abs(el(B,1,1)-el(GB,2,1))<1e-12
        test_assert "paired grid components ratio `method'/`distance'" abs(el(GB,4,2)-el(GY,4,2)/el(GD,4,2))<1e-12
        foreach sg in 4 8 {
            local row = cond(`sg'==4,2,4)
            clate_dnn py, denominator(pd) pseudooutcome xcols(x) ///
                x0matrix(CT) method(`method') distance(`distance') k(2) s(`sg')
            test_assert "scalar/grid estimate `method'/`distance'/`sg'" abs(_b[theta2]-el(GB,`row',2))<1e-12
            test_assert "scalar/grid SE `method'/`distance'/`sg'" abs(_se[theta2]-el(GS,`row',2))<1e-10
        }
    }
}
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0(0) tuningmode(grid)
test_assert "paper lambda default preserved" e(lambda)==.4
test_assert "paper grid defaults present" e(tuning_n_candidates)==4
