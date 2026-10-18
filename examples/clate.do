// Experimental CLATE with oracle instrument-arm predictions and paired ratios.
// Run after installing the package from master, or with src/
// on the adopath. This example illustrates mechanics, not fitted learners.
version 17
clear
set more off

// A deterministic paired-signal calculation, independent of causal claims.
// Common weights are (1/2,1/3,1/6,0). A deleted first stage is zero, but
// paired component covariance remains finite: no deleted ratio is formed.
set obs 4
gen double x = _n
gen double phi_y = 2 in 1
replace phi_y = 1 in 2
replace phi_y = 4 in 3
replace phi_y = 3 in 4
gen double phi_d = 1
replace phi_d = -1 in 2
replace phi_d = 2 in 3
clate_dnn phi_y, pseudooutcome denominator(phi_d) xcols(x) x0(0) s(2)
matrix fixture_components = e(components)
matrix fixture_covariance = e(component_vcov)
matrix fixture_ratio = e(theta_vcov)
assert abs(fixture_components[1,1]-2) < 1e-12
assert abs(fixture_components[1,2]-.5) < 1e-12
assert abs(_b[theta]-4) < 1e-12
assert abs(fixture_covariance[1,2]-.5) < 1e-12
assert abs(fixture_ratio[1,1]-34) < 1e-10

// An exactly zero aggregate first stage remains visible, never replaced.
gen double zero_stage = 0
clate_dnn phi_y, pseudooutcome denominator(zero_stage) ///
    xcols(x) x0(0) s(2)
matrix zero_components = e(components)
matrix zero_estimates = e(estimates)
matrix zero_status = e(status)
assert zero_components[1,2] == 0
assert missing(zero_estimates[1,1])
assert zero_status[1,1] == 1 & zero_status[1,5] == 0

// One-sided encouragement design with known instrument-arm predictions.
clear
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
clate_dnn Y, xcols(x1 x2) x0matrix(X0) dcol(D) acol(A) ecol(e) ///
    my0col(my0) my1col(my1) md1col(md1) onesided ///
    training(oracle) tuningmode(grid)
matrix oracle_estimates = e(estimates)
matrix oracle_covariance = e(theta_vcov)
matrix oracle_components = e(components)
matrix oracle_status = e(status)
matrix oracle_report = r(summary)
assert e(N) == 160 & e(m) == 3
assert rowsof(oracle_components) == 3 & colsof(oracle_components) == 2
assert rowsof(oracle_covariance) == 3 & colsof(oracle_covariance) == 3
forvalues j = 1/3 {
    assert oracle_status[`j',1] == 0 & oracle_status[`j',5] == 1
    assert abs(oracle_estimates[1,`j'] ///
        - oracle_components[`j',1]/oracle_components[`j',2]) < 1e-12
}
estimates store clate_oracle

// Reconstruct the paired signals with exactly the same propensity clipping.
gen double e_used = min(max(e,1e-8),1-1e-8)
gen double phi_y = my1-my0 + A*(Y-my1)/e_used ///
    - (1-A)*(Y-my0)/(1-e_used)
gen double phi_d = md1 + A*(D-md1)/e_used - (1-A)*D/(1-e_used)
clate_dnn phi_y, pseudooutcome denominator(phi_d) ///
    xcols(x1 x2) x0matrix(X0) training(oracle) tuningmode(grid)
matrix paired_estimates = e(estimates)
matrix paired_covariance = e(theta_vcov)
forvalues j = 1/3 {
    assert abs(paired_estimates[1,`j']-oracle_estimates[1,`j']) < 1e-12
    forvalues h = 1/3 {
        assert abs(paired_covariance[`j',`h'] ///
            - oracle_covariance[`j',`h']) < 1e-10
    }
}

// Restore the nuisance-based batch for instrument-overlap diagnostics/plots.
estimates restore clate_oracle
clate_dnn, level(90)
matrix oracle_report90 = r(summary)
matrix list oracle_components
matrix list oracle_status
clate_dnn_diagnostics, detail
matrix component_report = r(components)
assert rowsof(component_report) == 3 & colsof(component_report) == 9
clate_dnn_plot, type(estimates) name(clate_example_estimates) replace nodraw
clate_dnn_plot, type(overlap) view(mass) name(clate_example_overlap) replace nodraw
clate_dnn_plot, type(influence) name(clate_example_influence) replace nodraw
clate_dnn_plot, type(tuning) allcandidates name(clate_example_tuning) replace nodraw

// Leave the three-target CLATE fit active for interactive and install checks.
assert "`e(cmd)'" == "clate_dnn" & e(m) == 3
