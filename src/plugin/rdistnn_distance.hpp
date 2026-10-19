// Pure C++17 L1/L2 distances, derived from RdistNN/src/dist_l1.cpp and
// dist_l2.cpp with Rcpp removed. Shared by scalar, batch and legacy L2 callers.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace rdistnn {

enum class distance_metric { l1, l2 };

// Match Mata arithmetic selection before its finite missing-value boundary.
inline constexpr double distance_l2_range_limit=0x1p1023;

// Shared L1/L2 exceptional subtraction: ordinary values keep subtract-then-
// scale rounding; exact binary halves rescue an unrepresentable difference.
inline double scaled_distance_difference(double x,double target,double inverse)
{
    const double difference=x-target;
    if(std::fabs(difference)>=distance_l2_range_limit && std::fabs(inverse)<1) {
        const volatile double half_x=x*0.5,half_t=target*0.5;
        const double half_difference=half_x-half_t;
        return half_difference*(inverse+inverse);
    }
    return difference*inverse;
}

// R and Mata round the squared coordinate before adding it. An explicit
// store prevents FMA contraction from breaking stable ties on capable hosts.
inline double rounded_distance_square(double value)
{
    const volatile double product=value*value;
    return product;
}

// Revisit only exceptional L2 rows. The scale/ssq recurrence avoids squaring
// the largest coordinate; ordinary rows keep their original rounded arithmetic.
inline double distance_l2_range_row(const double* X, int n, int p, int i,
                                   const double* target, int stride_target,
                                   const double* inv_s, bool squared)
{
    double scale=0,ssq=1;
    for(int j=0;j<p;++j) {
        const double x=X[i+static_cast<std::size_t>(j)*n];
        const double t=target[static_cast<std::size_t>(j)*stride_target];
        double d=x-t;
        if(inv_s) d=scaled_distance_difference(x,t,inv_s[j]);
        const double a=std::fabs(d);
        if(!std::isfinite(a)) return a;
        if(a!=0) {
            if(scale<a) {
                const volatile double term=ssq*rounded_distance_square(scale/a);
                ssq=1+term;
                scale=a;
            } else ssq+=rounded_distance_square(a/scale);
        }
    }
    return squared?(scale*ssq)*scale:scale*std::sqrt(ssq);
}

// Validate dimensions/products before reading arrays or allocating a buffer.
// Bound both byte counts and pointer differences, including 32-bit hosts.
inline std::size_t checked_distance_size(int rows, int cols)
{
    if (rows<1 || cols<1)
        throw std::invalid_argument("distance matrix dimensions must be positive");
    const auto r=static_cast<std::size_t>(rows), c=static_cast<std::size_t>(cols);
    const auto max_elements=static_cast<std::size_t>(
        (std::numeric_limits<std::ptrdiff_t>::max)())/sizeof(double);
    if (r>max_elements/c)
        throw std::invalid_argument("distance matrix dimensions exceed addressable storage");
    return r*c;
}

inline void validate_distance_values(const double* values, std::size_t size)
{
    if (!values) throw std::invalid_argument("distance input must not be null");
    for (std::size_t i=0;i<size;++i)
        if (!std::isfinite(values[i]))
            throw std::invalid_argument("distance inputs must be finite");
}

// X: n x p, X0: m x p, output: n x m, all column-major. Optional inv_s
// contains p strictly positive inverse scales. m=1 supplies scalar distances.
// The coordinate accumulation order follows the R native implementations.
// Squared distances are supported only for the legacy L2 compatibility API.
// Invalid arguments and arithmetic overflow leave the output untouched.
inline void distance_batch_core(const double* X, int n, int p,
                                const double* X0, int m, const double* inv_s,
                                double* out, distance_metric metric,
                                bool squared=false)
{
    const auto nx=checked_distance_size(n,p), nt=checked_distance_size(m,p);
    const auto nd=checked_distance_size(n,m);
    if (metric!=distance_metric::l1 && metric!=distance_metric::l2)
        throw std::invalid_argument("distance metric must be L1 or L2");
    if (squared && metric!=distance_metric::l2)
        throw std::invalid_argument("squared distances require L2");
    if (!out) throw std::invalid_argument("distance output must not be null");
    validate_distance_values(X,nx);
    validate_distance_values(X0,nt);
    if (inv_s) {
        validate_distance_values(inv_s,static_cast<std::size_t>(p));
        for (int j=0;j<p;++j)
            if (inv_s[j]<=0)
                throw std::invalid_argument("inverse scales must be strictly positive");
    }
    // Stage values so a late overflow or an input/output alias is harmless.
    std::vector<double> result(nd);
    const auto stride_n=static_cast<std::size_t>(n), stride_m=static_cast<std::size_t>(m);
    for (int q=0;q<m;++q) for (int i=0;i<n;++i) {
        double acc=0;
        bool exceptional=false;
        for (int j=0;j<p;++j) {
            const double diff=X[i+static_cast<std::size_t>(j)*stride_n]-
                              X0[q+static_cast<std::size_t>(j)*stride_m];
            const double d=inv_s?scaled_distance_difference(
                X[i+static_cast<std::size_t>(j)*stride_n],
                X0[q+static_cast<std::size_t>(j)*stride_m],inv_s[j]):diff;
            if(metric==distance_metric::l1) acc+=std::fabs(d);
            else {
                const double square=rounded_distance_square(d);
                // Nonfinite squares propagate to acc; inspect it after the loop.
                exceptional=exceptional ||
                    (inv_s && std::fabs(diff)>=distance_l2_range_limit) ||
                    (d!=0 && square<std::numeric_limits<double>::min());
                acc+=square;
            }
        }
        const double value=metric==distance_metric::l2 ?
            ((exceptional || !std::isfinite(acc) || acc>=distance_l2_range_limit) ?
                distance_l2_range_row(X,n,p,i,X0+q,m,inv_s,squared) :
                (squared?acc:std::sqrt(acc))) : acc;
        if (!std::isfinite(value))
            throw std::invalid_argument("distance arithmetic overflowed");
        result[i+static_cast<std::size_t>(q)*stride_n]=value;
    }
    std::copy(result.begin(),result.end(),out);
}

} // namespace rdistnn
