// rdistnn_plugin.cpp
// Stata plugin entry point for RdistNN.
//
// Matrix operations: distance_api, distance_batch, and rank_loo.
// Current ado/Mata callers provide matrices and retain host-side rank tables.
//
// Compiling:
//   Include stplugin.c in the build so that _stata_ and pginit are defined.
//   See Makefile / Makefile.win in this directory.

#include "stplugin.h"
#include "rdistnn_core.hpp"

#include <climits>    // INT_MAX for rank_loo dimensions
#include <cmath>      // isfinite
#include <cstring>    // strcmp
#include <cstdio>     // snprintf
#include <vector>
#include <string>
#include <exception>

// Matrix API: plugin call rdistnn_plugin, rank_loo ranks signals table output
// Matrices are supplied by Mata; ranks n x m, signals n x q, table (n-1) x 1,
// preallocated output n x (m*q). No dataset or distance dependence.
static ST_retcode do_rank_loo(int argc, char *argv[])
{
    if (argc!=5) {
        SF_error((char *)"rdistnn_plugin rank_loo: require ranks signals table output matrix names.\n");
        return 198;
    }
    const int n=SF_row(argv[1]), m=SF_col(argv[1]), q=SF_col(argv[2]);
    if (n<2 || m<1 || q<1 || m>INT_MAX/q || SF_row(argv[2])!=n ||
        SF_row(argv[3])!=n-1 || SF_col(argv[3])!=1 ||
        SF_row(argv[4])!=n || SF_col(argv[4])!=m*q) {
        SF_error((char *)"rdistnn_plugin rank_loo: incompatible matrix dimensions.\n");
        return 198;
    }
    try {
        std::vector<int> ranks(static_cast<std::size_t>(n)*m);
        std::vector<double> signal(static_cast<std::size_t>(n)*q), table(n-1);
        ST_double value;
        ST_retcode rc;
        for (int j=0; j<m; ++j) for (int i=0; i<n; ++i) {
            if ((rc=SF_mat_el(argv[1],i+1,j+1,&value))) return rc;
            if (SF_is_missing(value) || !std::isfinite(value) || value<1 || value>n || value!=std::floor(value))
                throw std::invalid_argument("ranks must be integer permutations of 1..n");
            ranks[i+static_cast<std::size_t>(j)*n]=static_cast<int>(value);
        }
        for (int j=0; j<q; ++j) for (int i=0; i<n; ++i) {
            if ((rc=SF_mat_el(argv[2],i+1,j+1,&value))) return rc;
            if (SF_is_missing(value) || !std::isfinite(value))
                throw std::invalid_argument("signals must be finite and nonmissing");
            signal[i+static_cast<std::size_t>(j)*n]=value;
        }
        for (int i=0; i<n-1; ++i) {
            if ((rc=SF_mat_el(argv[3],i+1,1,&value))) return rc;
            if (SF_is_missing(value) || !std::isfinite(value))
                throw std::invalid_argument("table must be finite and nonmissing");
            table[i]=value;
        }
        std::vector<double> out(static_cast<std::size_t>(n)*m*q);
        rdistnn::rank_loo(ranks.data(),signal.data(),table.data(),n,m,q,out.data());
        for (double v: out) if (!std::isfinite(v) || SF_is_missing(v))
            throw std::invalid_argument("nonfinite prefix/suffix sum");
        for (int j=0; j<m*q; ++j) for (int i=0; i<n; ++i)
            if ((rc=SF_mat_store(argv[4],i+1,j+1,out[i+static_cast<std::size_t>(j)*n]))) return rc;
    }
    catch (const std::exception &error) {
        std::string msg="rdistnn_plugin rank_loo: "+std::string(error.what())+".\n";
        SF_error(const_cast<char *>(msg.c_str()));
        return 198;
    }
    return 0;
}

// Matrix API: distance_batch metric X X0 inverse_scales output.
// inverse_scales is a 1 x p matrix name or literal "none". Validate and read
// every input before publishing any cell, including when output aliases input.
static ST_retcode read_distance_matrix(char* name,int n,int p,std::vector<double>& values)
{
    values.resize(rdistnn::checked_distance_size(n,p));
    ST_double value;
    for (int j=0;j<p;++j) for (int i=0;i<n;++i) {
        const ST_retcode rc=SF_mat_el(name,i+1,j+1,&value);
        if (rc) return rc;
        if (SF_is_missing(value) || !std::isfinite(value))
            throw std::invalid_argument("distance inputs must be finite and nonmissing");
        values[i+static_cast<std::size_t>(j)*n]=value;
    }
    return 0;
}

static ST_retcode do_distance_batch(int argc,char* argv[])
{
    if (argc!=6)
        throw std::invalid_argument("distance_batch requires metric X X0 inverse_scales output");
    rdistnn::distance_metric metric;
    if (std::strcmp(argv[1],"l1")==0) metric=rdistnn::distance_metric::l1;
    else if (std::strcmp(argv[1],"l2")==0) metric=rdistnn::distance_metric::l2;
    else throw std::invalid_argument("distance_batch metric must be l1 or l2");
    const int n=SF_row(argv[2]), p=SF_col(argv[2]), m=SF_row(argv[3]);
    const auto nx=rdistnn::checked_distance_size(n,p);
    const auto nt=rdistnn::checked_distance_size(m,p);
    const auto nd=rdistnn::checked_distance_size(n,m);
    if (SF_col(argv[3])!=p || SF_row(argv[5])!=n || SF_col(argv[5])!=m)
        throw std::invalid_argument("distance_batch has incompatible matrix dimensions");
    const bool scaled=std::strcmp(argv[4],"none")!=0;
    if (scaled && (SF_row(argv[4])!=1 || SF_col(argv[4])!=p))
        throw std::invalid_argument("distance_batch inverse scales must be a 1 by p matrix");
    std::vector<double> X, X0, inverse;
    X.reserve(nx); X0.reserve(nt);
    ST_retcode rc=read_distance_matrix(argv[2],n,p,X);
    if (rc) return rc;
    rc=read_distance_matrix(argv[3],m,p,X0);
    if (rc) return rc;
    if (scaled) {
        rc=read_distance_matrix(argv[4],1,p,inverse);
        if (rc) return rc;
    }
    std::vector<double> out(nd);
    rdistnn::distance_batch_core(X.data(),n,p,X0.data(),m,
                                 scaled?inverse.data():nullptr,out.data(),metric);
    for (double v:out)
        if (SF_is_missing(v))
            throw std::invalid_argument("distance result exceeds Stata's finite range");
    for (int j=0;j<m;++j) for (int i=0;i<n;++i) {
        rc=SF_mat_store(argv[5],i+1,j+1,out[i+static_cast<std::size_t>(j)*n]);
        if (rc) return rc;
    }
    return 0;
}

// Capability probe for the shared distance dispatcher. A compatible plugin
// writes API version 1 to a preallocated 1 x 1 matrix.
static ST_retcode do_distance_api(int argc,char* argv[])
{
    if (argc!=2 || SF_row(argv[1])!=1 || SF_col(argv[1])!=1)
        throw std::invalid_argument("distance_api requires a preallocated 1 by 1 output matrix");
    return SF_mat_store(argv[1],1,1,1.0);
}

// ── Plugin entry point ─────────────────────────────────────────────────────

STDLL stata_call(int argc, char *argv[])
{
    // Includes allocation failures in all handlers: no C++ exception may cross
    // the Stata plugin ABI. Fixed buffers keep error handling allocation-free.
    try {
        if (argc<1 || !argv)
            throw std::invalid_argument("function name required as first argument");
        for (int i=0;i<argc;++i)
            if (!argv[i] || argv[i][0]=='\0')
                throw std::invalid_argument("plugin arguments must be nonempty");
        if (std::strcmp(argv[0],"rank_loo")==0) return do_rank_loo(argc,argv);
        if (std::strcmp(argv[0],"distance_batch")==0) return do_distance_batch(argc,argv);
        if (std::strcmp(argv[0],"distance_api")==0) return do_distance_api(argc,argv);
        throw std::invalid_argument("unknown operation; supported: rank_loo distance_batch distance_api");
    } catch (const std::exception& error) {
        char msg[768];
        std::snprintf(msg,sizeof(msg),"rdistnn_plugin: %s.\n",error.what());
        SF_error(msg);
        return 198;
    } catch (...) {
        SF_error((char*)"rdistnn_plugin: unexpected native failure.\n");
        return 198;
    }
}
