// dev-setup.do
// Run once per machine to add the src/ directory to Stata's adopath so
// that cate_dnn, np_dnn, weights_dnn etc. are found during development.
//
// Usage (from Stata command line or do-file editor):
//   do "D:\Jakob_Clouds\Dropbox\Research\Individual_Projects\U_Statistics\02_Code\RdistNN_stata\dev-setup.do"
//
// To make permanent across sessions, paste the adopath line into your
// profile.do (found via: di c(sysdir_personal)).

local pkg_src `"D:\Jakob_Clouds\Dropbox\Research\Individual_Projects\U_Statistics\02_Code\RdistNN_stata\src"'

// Check src/ exists before modifying adopath
capture confirm file `"`pkg_src'/_rdistnn_mata.ado"'
if _rc {
    di as error "src/ directory not found at expected path."
    di as error "Edit the local pkg_src path in this script."
    exit 601
}

adopath ++ `"`pkg_src'"'
di as result "rdistnn src/ added to adopath for this session."
di as result "Current adopath:"
adopath
