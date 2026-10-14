// dev-setup.do
// Run once per session to add the src/ directory to Stata's adopath so
// that cate_dnn, np_dnn, weights_dnn etc. are found during development.
//
// Usage: cd to the repo root, then:
//   do dev-setup.do
//
// To make permanent across sessions, paste the adopath line into your
// profile.do (found via: di c(sysdir_personal)).

// Resolve path relative to the current working directory so this script
// works on any machine without hardcoded paths.
// Forward slashes also keep Windows UNC/WSL paths usable by ado discovery.
local here = subinstr("`c(pwd)'", "\", "/", .)
local pkg_src `"`here'/src"'

// Check src/ exists before modifying adopath
capture confirm file `"`pkg_src'/_rdistnn_mata.ado"'
if _rc {
    di as error "src/ directory not found at `pkg_src'."
    di as error "Run this script from the repository root (cd to the folder"
    di as error "that contains dev-setup.do and the src/ subdirectory)."
    exit 601
}

adopath ++ `"`pkg_src'"'
di as result "rdistnn src/ added to adopath for this session."
di as result "Current adopath:"
adopath
