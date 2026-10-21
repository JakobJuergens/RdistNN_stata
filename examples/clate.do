// Conditional LATE with oracle instrument-arm predictions.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This self-contained example clears the active dataset.
version 17
set more off
clear
// This one-sided encouragement DGP illustrates inputs, not fitted learners.
set obs 160
set seed 4711
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double e = invlogit(.3*x1-.2*x2)
gen byte A = runiform()<e
gen double my0 = .5*x1-.2*x2
gen double md1 = invlogit(.4+.2*x1)
gen double md0 = 0
gen double tau = 1+.3*sin(x1)
gen byte D = A*(runiform()<md1)
gen double Y = my0 + tau*D + rnormal()*.3
gen double my1 = my0 + md1*tau
matrix X0 = (-.5,0 \ 0,0 \ .5,0)

// Arm means condition on instrument A. The command does not fit nuisances.
clate_dnn Y, xcols(x1 x2) x0matrix(X0) dcol(D) acol(A) ecol(e) ///
    my0col(my0) my1col(my1) md1col(md1) onesided ///
    training(oracle) tuningmode(grid)
matrix components = e(components)
matrix statuses = e(status)
matrix effect_report = r(summary)
estimates store clate_oracle
matrix list components
matrix list effect_report
matrix list statuses
clate_dnn_diagnostics

// Supplied paired signals produce the same effects on the same sample.
gen double e_used = min(max(e,1e-8),1-1e-8)
gen double phi_y = my1-my0 + A*(Y-my1)/e_used ///
    - (1-A)*(Y-my0)/(1-e_used)
gen double phi_d = md1 + A*(D-md1)/e_used - (1-A)*D/(1-e_used)
clate_dnn phi_y, pseudooutcome denominator(phi_d) ///
    xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix paired_report = r(summary)
matrix list paired_report

// Restore the nuisance-based fit for instrument-overlap diagnostics.
// Inspect the first stage: finite Wald bounds are not a weak-IV guarantee.
estimates restore clate_oracle
clate_dnn, level(90)
clate_dnn_plot, type(estimates) name(clate_example_estimates) replace nodraw
clate_dnn_plot, type(overlap) view(mass) name(clate_example_overlap) replace nodraw
clate_dnn_plot, type(influence) targets(2) name(clate_example_influence) replace nodraw
clate_dnn_plot, type(tuning) allcandidates name(clate_example_tuning) replace nodraw
