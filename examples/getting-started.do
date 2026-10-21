// A first CATE fit and a conditional mean, using oracle nuisances.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This self-contained example clears the active dataset.
version 17
set more off
clear
set obs 300
set seed 42
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double e_oracle = invlogit(0.5*x1 - 0.5*x2)
gen byte D = runiform() < e_oracle
gen double m0_oracle = 1 + x1 - x2
gen double m1_oracle = m0_oracle + 0.5 + 0.5*sin(x1)
gen double Y = m0_oracle + D*(m1_oracle-m0_oracle) + rnormal()*0.5

// In applications, replace these oracle columns with aligned external predictions.
matrix X0 = (0,0 \ 0.5,0 \ -0.5,0)
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle)
matrix cate_report = r(summary)
estimates store first_cate
matrix list cate_report
cate_dnn_diagnostics
cate_dnn_plot, type(estimates) name(cate_example) replace nodraw

// NP estimates E[Y|X], using the same targets without nuisance inputs.
np_dnn Y, xcols(x1 x2) x0matrix(X0)
matrix mean_report = r(summary)
matrix list mean_report
np_dnn_diagnostics

// Restore the matching fit before CATE postestimation.
estimates restore first_cate
cate_dnn, level(90)
