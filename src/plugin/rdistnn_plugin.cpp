// rdistnn_plugin.cpp
// Stata plugin entry point for RdistNN.
//
// Plugin call convention (from Stata .ado):
//
//   plugin call rdistnn_plugin [kx_vars] [w_out_var] [if] [in] , ///
//       weights_l2 [method] [k] [s] [x0_1] ... [x0_kx]
//
// Variables passed to plugin call (1-based indices via SF_vdata/SF_vstore):
//   1 .. kx   : covariate columns (read)
//   kx + 1    : output weight column (written)
//
// argv[] after the comma (0-based):
//   argv[0]           : function name, currently only "weights_l2"
//   argv[1]           : method: "exact" or "asymp"
//   argv[2]           : k  (integer, # selected within each s-subset)
//   argv[3]           : s  (integer, subset/kernel order >= k; larger is tighter)
//   argv[4 .. 4+kx-1] : x0 components (one per covariate)
//
// Compiling:
//   Include stplugin.c in the build so that _stata_ and pginit are defined.
//   See Makefile / Makefile.win in this directory.

#include "stplugin.h"
#include "rdistnn_core.hpp"

#include <cerrno>     // errno, ERANGE
#include <climits>    // INT_MIN, INT_MAX
#include <cmath>      // isfinite
#include <cstdlib>    // strtol, strtod
#include <cstring>    // strcmp
#include <cstdio>     // snprintf
#include <vector>
#include <string>
#include <exception>

// ── helpers ────────────────────────────────────────────────────────────────

// Build the list of in-sample (if/in) observation numbers (1-based).
static std::vector<ST_int> collect_obs()
{
    std::vector<ST_int> obs;
    obs.reserve(static_cast<std::size_t>(SF_nobs()));
    ST_int lo = SF_in1(), hi = SF_in2();
    for (ST_int o = lo; o <= hi; ++o) {
        if (SF_ifobs(o)) obs.push_back(o);
    }
    return obs;
}

// Read column j (1-based Stata variable index) for all obs_list entries
// into out[0..n-1].
static ST_retcode read_column(int j_stata, const std::vector<ST_int> &obs_list,
                               double *out)
{
    for (std::size_t idx = 0; idx < obs_list.size(); ++idx) {
        ST_double z;
        ST_retcode rc = SF_vdata(j_stata, obs_list[idx], &z);
        if (rc) return rc;
        if (SF_is_missing(z) || !std::isfinite(static_cast<double>(z))) {
            SF_error((char *)"rdistnn_plugin weights_l2: covariates must contain only finite, nonmissing values in the marked sample.\n");
            return 198;
        }
        out[idx] = static_cast<double>(z);
    }
    return 0;
}

static bool parse_int_arg(const char *arg, const char *name, int *out)
{
    if (arg == nullptr || arg[0] == '\0') {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "rdistnn_plugin weights_l2: %s must be a nonempty integer.\n", name);
        SF_error(msg);
        return false;
    }

    errno = 0;
    char *end = nullptr;
    const long value = std::strtol(arg, &end, 10);
    if (end == arg || *end != '\0' || errno == ERANGE ||
        value < INT_MIN || value > INT_MAX) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "rdistnn_plugin weights_l2: %s must be an integer, got '%s'.\n",
            name, arg);
        SF_error(msg);
        return false;
    }

    *out = static_cast<int>(value);
    return true;
}

static bool parse_double_arg(const char *arg, const char *name, double *out)
{
    if (arg == nullptr || arg[0] == '\0') {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "rdistnn_plugin weights_l2: %s must be a nonempty number.\n", name);
        SF_error(msg);
        return false;
    }

    errno = 0;
    char *end = nullptr;
    const double value = std::strtod(arg, &end);
    if (end == arg || *end != '\0' || errno == ERANGE || !std::isfinite(value)) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "rdistnn_plugin weights_l2: %s must be a finite number, got '%s'.\n",
            name, arg);
        SF_error(msg);
        return false;
    }

    *out = value;
    return true;
}

// ── "weights_l2" handler ───────────────────────────────────────────────────
//
// Reads kx covariate columns + x0 from argv, computes DNN L2 weights,
// and writes them to the last Stata variable in the varlist.

static ST_retcode do_weights_l2(int argc, char *argv[])
{
    // Number of Stata variables in the plugin call varlist
    const int nvars = static_cast<int>(SF_nvars());
    if (nvars < 2) {
        SF_error((char *)"rdistnn_plugin weights_l2: kx must be >= 1 (provide at least one covariate and one output variable).\n");
        return 198;
    }
    const int kx   = nvars - 1;   // # covariates
    const int wvar = nvars;        // 1-based index of the output weight variable

    // argc must be exactly: function_name, method, k, s, x0_1, ..., x0_kx
    const int expected_argc = 4 + kx;
    if (argc != expected_argc) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "rdistnn_plugin weights_l2: need exactly %d arguments (weights_l2 method k s x0_1 ... x0_%d), got %d.\n",
            expected_argc, kx, argc);
        SF_error(msg);
        return 198;
    }

    // argv[1]: method ("exact" or "asymp")
    const bool asymp = (std::strcmp(argv[1], "asymp") == 0);
    if (!asymp && std::strcmp(argv[1], "exact") != 0) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "rdistnn_plugin weights_l2: unknown method '%s'. Use 'exact' or 'asymp'.\n", argv[1]);
        SF_error(msg);
        return 198;
    }

    int k = 0;
    int s = 0;
    if (!parse_int_arg(argv[2], "k", &k)) return 198;
    if (!parse_int_arg(argv[3], "s", &s)) return 198;

    std::vector<double> x0(kx);
    for (int j = 0; j < kx; ++j) {
        char name[32];
        std::snprintf(name, sizeof(name), "x0_%d", j + 1);
        if (!parse_double_arg(argv[4 + j], name, &x0[j])) return 198;
    }

    // Collect in-sample observations (respects if/in + touse variable)
    const std::vector<ST_int> obs_list = collect_obs();
    if (obs_list.size() > static_cast<std::size_t>(INT_MAX)) {
        SF_error((char *)"rdistnn_plugin weights_l2: marked sample is too large for the plugin integer interface.\n");
        return 198;
    }
    const int n = static_cast<int>(obs_list.size());
    if (n < 1) {
        SF_error((char *)"rdistnn_plugin weights_l2: n must be >= 1; the marked sample is empty.\n");
        return 198;
    }
    if (kx < 1) {
        SF_error((char *)"rdistnn_plugin weights_l2: kx must be >= 1.\n");
        return 198;
    }
    if (k < 1) {
        SF_error((char *)"rdistnn_plugin weights_l2: k must be >= 1.\n");
        return 198;
    }
    if (s < k) {
        SF_error((char *)"rdistnn_plugin weights_l2: s must satisfy s >= k.\n");
        return 198;
    }
    if (s > n) {
        SF_error((char *)"rdistnn_plugin weights_l2: s must satisfy s <= n in the marked sample.\n");
        return 198;
    }

    // Read X: n x kx, column-major (column j stored at X + j*n)
    std::vector<double> X(static_cast<std::size_t>(n) * kx);
    for (int j = 0; j < kx; ++j) {
        ST_retcode rc = read_column(j + 1, obs_list, X.data() + static_cast<std::size_t>(j) * n);
        if (rc) return rc;
    }

    // Compute DNN weights via the L2 pipeline
    std::vector<double> w(n);
    try {
        rdistnn::dnn_weights_l2(X.data(), n, kx, x0.data(),
                                 /*inv_s=*/nullptr, k, s, w.data(), asymp);
    }
    catch (const std::exception &error) {
        char msg[512];
        std::snprintf(msg, sizeof(msg),
            "rdistnn_plugin weights_l2: invalid core arguments: %s.\n",
            error.what());
        SF_error(msg);
        return 198;
    }

    // Write weights back to the output variable in the Stata dataset
    for (int i = 0; i < n; ++i) {
        ST_retcode rc = SF_vstore(wvar, obs_list[i], static_cast<ST_double>(w[i]));
        if (rc) return rc;
    }

    return 0;
}

// ── Plugin entry point ─────────────────────────────────────────────────────

STDLL stata_call(int argc, char *argv[])
{
    if (argc < 1) {
        SF_error((char *)"rdistnn_plugin: function name required as first argument.\n");
        return 198;
    }

    if (std::strcmp(argv[0], "weights_l2") == 0) {
        return do_weights_l2(argc, argv);
    }

    char msg[256];
    std::snprintf(msg, sizeof(msg),
        "rdistnn_plugin: unknown function '%s'. Supported: weights_l2\n", argv[0]);
    SF_error(msg);
    return 198;
}
