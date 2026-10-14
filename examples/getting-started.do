// RdistNN Stata onboarding: simulated data and explicitly oracle nuisances.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This example clears the active dataset.
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

matrix X0 = (0,0 \ 0.5,0 \ -0.5,0)
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid)
matrix list e(b)
matrix list e(V)
cate_dnn_diagnostics
cate_dnn_plot, type(estimates) name(cate_example) replace nodraw
cate_dnn_plot, type(tuning) allcandidates name(cate_tuning_example) replace nodraw
cate_dnn_plot, type(overlap) view(mass) name(cate_overlap_example) replace nodraw

np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
np_dnn_diagnostics
np_dnn_plot, type(estimates) name(np_example) replace nodraw
np_dnn_plot, type(outcome) name(np_outcome_example) replace nodraw

// Precomputed oracle AIPW score: estimate-only CATE at one target.
gen double eclip = min(max(e_oracle,1e-8),1-1e-8)
gen double rho = m1_oracle-m0_oracle + D*(Y-m1_oracle)/eclip ///
    - (1-D)*(Y-m0_oracle)/(1-eclip)
cate_dnn rho, xcols(x1 x2) x0(0 0) pseudooutcome nose

// Robinson/R-learner route with the oracle pooled outcome mean.
gen double m_oracle = e_oracle*m1_oracle + (1-e_oracle)*m0_oracle
cate_dnn Y, tcol(D) xcols(x1 x2) x0(0 0) ecol(e_oracle) mcol(m_oracle)

weights_dnn x1 x2, x0(0 0) k(1) generate(w_example)
summarize w_example
