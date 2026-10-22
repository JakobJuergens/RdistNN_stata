version 17
capture mata: mata drop rdistnn_test_finite_close()
mata:
real scalar rdistnn_test_finite_close(real matrix actual, real matrix expected,
    real scalar tolerance)
{
    real matrix difference
    if (rows(actual)!=rows(expected) | cols(actual)!=cols(expected) |
        rows(actual)==0 | cols(actual)==0) return(0)
    if (hasmissing(actual) | hasmissing(expected)) return(0)
    difference=actual-expected
    if (hasmissing(difference)) return(0)
    if (tolerance==0) return(max(abs(difference))==0)
    return(max(abs(difference))<tolerance)
}
end
