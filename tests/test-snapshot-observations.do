version 17
do "tests/helpers/fitted-data.do"
test_section "snapshot plots and numeric contribution extraction"
foreach family in np cate clate {
    local opts "xcols(fr_x) x0matrix(FR_NEW) s(4)"
    if "`family'"=="np" quietly np_dnn fr_y, `opts'
    if "`family'"=="cate" quietly cate_dnn fr_y, pseudooutcome `opts'
    if "`family'"=="clate" quietly clate_dnn fr_y, pseudooutcome denominator(fr_den) `opts'
    quietly `family'_dnn_plot, type(weights) view(kernel) targets(3 1) name(fr_kernel) replace nodraw
    matrix FR_KERNEL=r(plot_data)
    if "`family'"!="np" {
        quietly rdistnn_results, type(influence) targets(3 1)
        matrix FR_INFLUENCE=r(data)
    }
    preserve
    clear
    quietly `family'_dnn_plot, type(weights) view(kernel) targets(3 1) name(fr_kernel) replace nodraw
    matrix FR_COLD=r(plot_data)
    mata: st_numscalar("_fr_error",max(abs(st_matrix("FR_KERNEL")-st_matrix("FR_COLD"))))
    test_assert "`family': kernel graph uses saved geometry after clear" scalar(_fr_error)==0 & _N==0
    if "`family'"!="np" {
        quietly `family'_dnn_plot, type(influence) targets(3 1) name(fr_influence) replace nodraw
        matrix FR_COLD=r(plot_data)
        mata: st_numscalar("_fr_error",max(abs(st_matrix("FR_INFLUENCE")-st_matrix("FR_COLD"))))
        test_assert "`family': influence graph matches graph-free numeric extraction after clear" scalar(_fr_error)==0 & _N==0
    }
    restore
}
