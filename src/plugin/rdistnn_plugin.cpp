// rdistnn_plugin.cpp
// Stata plugin entry point for RdistNN.
//
// Matrix operations: distance_api, distance_batch, loo_packed_api and loo_packed.
// Current ado/Mata callers provide matrices and retain host-side rank tables.
//
// Compiling:
//   Include stplugin.c in the build so that _stata_ and pginit are defined.
//   See Makefile / Makefile.win in this directory.

#include "stplugin.h"
#include "rdistnn_core.hpp"

#include <climits>    // INT_MAX for validated payload dimensions
#include <cmath>      // isfinite
#include <cstring>    // strcmp
#include <cstdio>     // snprintf
#include <vector>
#include <string>
#include <exception>

// Near-square staging prevents default Stata row/column stripe construction
// from dominating long, thin matrices. Physical cells are row-major; payload
// arrays retain the existing logical column-major kernel layout.
struct packed_shape {
    int rows,cols;
    std::size_t capacity;
};

static std::size_t packed_count(std::size_t a,std::size_t b)
{
    const std::size_t exact_limit=static_cast<std::size_t>(9007199254740991ULL);
    if (b && a>std::numeric_limits<std::size_t>::max()/b)
        throw std::invalid_argument("packed matrix size overflows native dimensions");
    const std::size_t count=a*b;
    if (count>exact_limit)
        throw std::invalid_argument("packed matrix size exceeds exact host dimensions");
    return count;
}

static std::size_t packed_add(std::size_t a,std::size_t b)
{
    const std::size_t exact_limit=static_cast<std::size_t>(9007199254740991ULL);
    if (a>exact_limit || b>exact_limit-a)
        throw std::invalid_argument("packed matrix size exceeds exact host dimensions");
    return a+b;
}

static packed_shape packed_dimensions(char* name)
{
    const int rows=SF_row(name),cols=SF_col(name);
    if (rows<1 || cols<1 || rows<cols || rows-cols>1)
        throw std::invalid_argument("packed matrices must have balanced positive dimensions");
    return {rows,cols,packed_count(static_cast<std::size_t>(rows),cols)};
}

static void validate_packed_capacity(const packed_shape& shape,std::size_t length)
{
    if (shape.capacity<length || shape.capacity-length>=static_cast<std::size_t>(shape.rows))
        throw std::invalid_argument("packed matrix has incompatible payload dimensions");
}

static int packed_integer(double value,const char* name,int minimum)
{
    if (SF_is_missing(value) || !std::isfinite(value) || value<minimum ||
        value>INT_MAX || value!=std::floor(value))
        throw std::invalid_argument(std::string(name)+" must be supported integer dimensions");
    return static_cast<int>(value);
}

static ST_retcode do_loo_packed(int argc,char* argv[])
{
    if (argc!=3)
        throw std::invalid_argument("loo_packed requires input and output matrix names");
    const packed_shape input=packed_dimensions(argv[1]);
    if (input.capacity<5)
        throw std::invalid_argument("packed deletion metadata is incomplete");
    auto read=[&](std::size_t k,double& value) {
        return SF_mat_el(argv[1],static_cast<int>(k/input.cols)+1,
                         static_cast<int>(k%input.cols)+1,&value);
    };
    double header[5];ST_retcode rc;
    for (std::size_t k=0;k<5;++k) if ((rc=read(k,header[k]))) return rc;
    if (header[0]!=1 || (header[4]!=0 && header[4]!=1))
        throw std::invalid_argument("unsupported packed deletion version or group flag");
    const int n=packed_integer(header[1],"n",2),m=packed_integer(header[2],"m",1),
              q=packed_integer(header[3],"q",1);
    const bool grouped=header[4]==1;
    if (m>INT_MAX/q)
        throw std::invalid_argument("packed deletion target/component dimensions overflow");
    const std::size_t nr=packed_count(static_cast<std::size_t>(n),m),
                      ns=packed_count(static_cast<std::size_t>(n),q),
                      no=packed_count(nr,q);
    std::size_t length=packed_add(packed_add(packed_add(5,nr),ns),n-1);
    if (grouped) length=packed_add(length,n);
    validate_packed_capacity(input,length);
    const packed_shape output=packed_dimensions(argv[2]);
    validate_packed_capacity(output,no);
    // All dimensions and capacities are checked before native allocations.
    std::vector<int> ranks(nr),groups(grouped?n:0);
    std::vector<double> signal(ns),table(n-1),out(no);
    std::size_t cursor=5;double value;
    for (std::size_t k=0;k<nr;++k) {
        if ((rc=read(cursor++,value))) return rc;
        ranks[k]=packed_integer(value,"ranks",1);
        if (ranks[k]>n) throw std::invalid_argument("ranks must be permutations of 1..n");
    }
    for (double& v:signal) {
        if ((rc=read(cursor++,v))) return rc;
        if (SF_is_missing(v) || !std::isfinite(v))
            throw std::invalid_argument("packed signals must be finite and nonmissing");
    }
    for (double& v:table) {
        if ((rc=read(cursor++,v))) return rc;
        if (SF_is_missing(v) || !std::isfinite(v))
            throw std::invalid_argument("packed weights must be finite and nonmissing");
    }
    for (int& group:groups) {
        if ((rc=read(cursor++,value))) return rc;
        group=packed_integer(value,"groups",1);
        if (group>n) throw std::invalid_argument("groups must be identifiers in 1..n");
    }
    for (;cursor<input.capacity;++cursor) {
        if ((rc=read(cursor,value))) return rc;
        if (value!=0) throw std::invalid_argument("packed input padding must be zero");
    }
    if (grouped) rdistnn::block_rank_loo(ranks.data(),signal.data(),table.data(),groups.data(),n,m,q,out.data());
    else rdistnn::rank_loo(ranks.data(),signal.data(),table.data(),n,m,q,out.data());
    for (double v:out)
        if (SF_is_missing(v) || !std::isfinite(v))
            throw std::overflow_error("packed deletion result is outside the numeric range");
    // Materializing every input first also makes conformable aliases safe.
    for (std::size_t k=0;k<output.capacity;++k) {
        rc=SF_mat_store(argv[2],static_cast<int>(k/output.cols)+1,
                         static_cast<int>(k%output.cols)+1,k<no?out[k]:0.0);
        if (rc) return rc;
    }
    return 0;
}

static ST_retcode do_loo_packed_api(int argc,char* argv[])
{
    if (argc!=2 || SF_row(argv[1])!=1 || SF_col(argv[1])!=1)
        throw std::invalid_argument("loo_packed_api requires a 1 by 1 output matrix");
    return SF_mat_store(argv[1],1,1,1.0);
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
            throw std::overflow_error("distance result exceeds Stata's finite range");
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
        if (std::strcmp(argv[0],"loo_packed")==0) return do_loo_packed(argc,argv);
        if (std::strcmp(argv[0],"loo_packed_api")==0) return do_loo_packed_api(argc,argv);
        if (std::strcmp(argv[0],"distance_batch")==0) return do_distance_batch(argc,argv);
        if (std::strcmp(argv[0],"distance_api")==0) return do_distance_api(argc,argv);
        throw std::invalid_argument("unknown operation; supported: loo_packed loo_packed_api distance_batch distance_api");
    } catch (const std::bad_alloc&) {
        SF_error((char*)"rdistnn_plugin: insufficient memory; reduce the batch size.\n");
        return 909;
    } catch (const std::overflow_error&) {
        SF_error((char*)"rdistnn_plugin: numerical range exceeded; rescale covariates or signals.\n");
        return 498;
    } catch (const std::invalid_argument& error) {
        char msg[768];
        std::snprintf(msg,sizeof(msg),"rdistnn_plugin: %s.\n",error.what());
        SF_error(msg);
        return 198;
    } catch (const std::exception&) {
        SF_error((char*)"rdistnn_plugin: unexpected native execution failure; report a reproducible example.\n");
        return 498;
    } catch (...) {
        SF_error((char*)"rdistnn_plugin: unexpected native failure; report a reproducible example.\n");
        return 498;
    }
}
