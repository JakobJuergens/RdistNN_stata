// rdistnn_core.hpp
// Pure C++17 computation routines for the RdistNN Stata plugin.
// No Rcpp, no R dependencies.
//
// Provides:
//   rdistnn::dist_l2_core()       — row-wise L2 distances (column-major X)
//   rdistnn::argsort_indices()    — stable argsort (radix fast path)
//   rdistnn::ranks_from_values()  — 1-based ranks
//   rdistnn::dnn_weights_from_ranks() — DNN weight formula
//   rdistnn::dnn_weights_l2()     — full pipeline: dist -> ranks -> weights
//   rdistnn::rank_loo()           — prefix/suffix frozen-rank deletions
//
// Ported from RdistNN/src/{dist_l2.cpp, fast_order.cpp} with Rcpp removed.

#pragma once

#include <vector>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

namespace rdistnn {

inline void validate_sample_dimensions(int n, int p)
{
    if (n < 1) throw std::invalid_argument("n must be >= 1");
    if (p < 1) throw std::invalid_argument("p must be >= 1");
}

inline void validate_weight_parameters(int n, int k, int s)
{
    if (n < 1) throw std::invalid_argument("n must be >= 1");
    if (k < 1) throw std::invalid_argument("k must be >= 1");
    if (s < k) throw std::invalid_argument("s must be >= k");
    if (s > n) throw std::invalid_argument("s must be <= n");
}

inline void validate_ranks(const int *ranks, int n)
{
    if (ranks == nullptr) throw std::invalid_argument("ranks must not be null");
    std::vector<unsigned char> seen(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        const int rank = ranks[i];
        if (rank < 1 || rank > n || seen[static_cast<std::size_t>(rank - 1)] != 0)
            throw std::invalid_argument("ranks must be a permutation of 1..n");
        seen[static_cast<std::size_t>(rank - 1)] = 1;
    }
}

inline void validate_finite_array(const double *x, std::size_t size,
                                  const char *name)
{
    if (x == nullptr) throw std::invalid_argument(std::string(name) + " must not be null");
    for (std::size_t i = 0; i < size; ++i) {
        if (!std::isfinite(x[i]))
            throw std::invalid_argument(std::string(name) + " must contain only finite values");
    }
}

// ── L2 distance computation ────────────────────────────────────────────────
//
// X is n x p stored column-major (Stata and R both use column-major).
// x0 is length p.  inv_s[j] = 1/scale[j] for scaled distance; pass nullptr
// for unscaled Euclidean.  Writes n distances into out[].

inline void dist_l2_core(const double *X, int n, int p,
                          const double *x0,
                          const double *inv_s,   // nullptr = unscaled
                          double *out,
                          bool squared = false)
{
    validate_sample_dimensions(n, p);
    validate_finite_array(X, static_cast<std::size_t>(n) * p, "X");
    validate_finite_array(x0, static_cast<std::size_t>(p), "x0");
    if (inv_s != nullptr) validate_finite_array(inv_s, static_cast<std::size_t>(p), "inv_s");
    if (out == nullptr) throw std::invalid_argument("out must not be null");
    for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (int j = 0; j < p; ++j) {
            const double diff = X[i + (std::size_t)j * n] - x0[j];
            const double d    = inv_s ? diff * inv_s[j] : diff;
            acc += d * d;
        }
        out[i] = squared ? acc : std::sqrt(acc);
    }
}

// ── Stable argsort — radix fast path for non-negative finite values ────────

inline bool all_finite_nonneg(const double *x, int n) noexcept
{
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(x[i]) || x[i] < 0.0) return false;
    }
    return true;
}

namespace detail {

inline void argsort_radix64_nonneg(const double *x, int n, int *out)
{
    if (n <= 1) { if (n == 1) out[0] = 0; return; }

    std::vector<std::uint64_t> keys(n);
    std::memcpy(keys.data(), x, sizeof(double) * (std::size_t)n);

    for (int i = 0; i < n; ++i) out[i] = i;
    std::vector<int> buf(n);
    int *from = out, *to = buf.data();

    static constexpr int PASS_BITS = 11;
    static constexpr int RAD       = 1 << PASS_BITS;   // 2048
    static constexpr int PASSES    = (64 + PASS_BITS - 1) / PASS_BITS; // 6
    const std::uint64_t  MASK      = (std::uint64_t)(RAD - 1);
    std::vector<std::size_t> counts(RAD);

    for (int pass = 0; pass < PASSES; ++pass) {
        std::fill(counts.begin(), counts.end(), 0);
        const int shift = pass * PASS_BITS;
        for (int i = 0; i < n; ++i)
            ++counts[(keys[from[i]] >> shift) & MASK];
        std::size_t sum = 0;
        for (int b = 0; b < RAD; ++b) {
            const std::size_t c = counts[b]; counts[b] = sum; sum += c;
        }
        for (int i = 0; i < n; ++i) {
            const unsigned b = (unsigned)((keys[from[i]] >> shift) & MASK);
            to[counts[b]++] = from[i];
        }
        std::swap(from, to);
    }
    if (from != out)
        std::memcpy(out, from, sizeof(int) * (std::size_t)n);
}

struct StableCmp {
    const double *x;
    explicit StableCmp(const double *px) : x(px) {}
    static int nf_rank(double v) {
        return std::isnan(v) ? 2 : (v > 0.0 ? 1 : 0);
    }
    bool operator()(int i, int j) const {
        const double a = x[i], b = x[j];
        const bool ai = std::isfinite(a), aj = std::isfinite(b);
        if (ai && aj) { if (a < b) return true; if (b < a) return false; return i < j; }
        if (ai != aj)  return ai;
        const int ra = nf_rank(a), rb = nf_rank(b);
        return ra != rb ? ra < rb : i < j;
    }
};

} // namespace detail

// Stable argsort: out[k] = 0-based index of the k-th smallest value.
inline void argsort_indices(const double *x, int n, int *out)
{
    if (n < 1) throw std::invalid_argument("n must be >= 1");
    if (x == nullptr) throw std::invalid_argument("x must not be null");
    if (out == nullptr) throw std::invalid_argument("out must not be null");
    if (all_finite_nonneg(x, n)) {
        detail::argsort_radix64_nonneg(x, n, out);
    } else {
        for (int i = 0; i < n; ++i) out[i] = i;
        std::stable_sort(out, out + n, detail::StableCmp(x));
    }
}

// 1-based ranks: ranks_out[i] = rank of x[i] in ascending order (1 = smallest).
// Ties broken by original index.
inline void ranks_from_values(const double *x, int n, int *ranks_out)
{
    if (n < 1) throw std::invalid_argument("n must be >= 1");
    if (x == nullptr) throw std::invalid_argument("x must not be null");
    if (ranks_out == nullptr) throw std::invalid_argument("ranks_out must not be null");
    std::vector<int> ord(n);
    argsort_indices(x, n, ord.data());
    for (int pos = 0; pos < n; ++pos) ranks_out[ord[pos]] = pos + 1;
}

// ── DNN weight computation (Steele closed form / full WOR hypergeometric) ───
//
// For obs with global rank r (1 = nearest):
//
//   w_r = (1/k) * Σ_{l=1}^{k} C(r-1, l-1) · C(n-r, s-l) / C(n, s)
//
// Positive-weight support: r = 1, …, n-s+k  (NOT 1..s).  Larger s therefore
// shrinks support and tightens localization; for k=1 the edge is n-s+1.
// Here k counts the observations selected within each s-subset, not the global
// support size.  For k=1, w_r = C(n-r, s-1) / C(n, s) (Steele formula).
//
// Derivation: w_r is the marginal probability that the obs at global rank r
// is among the k nearest in a uniformly-random s-subset of all n obs.

// log C(n, k) via lgamma for numerical stability with large n.
inline double log_binom(int nn, int kk)
{
    if (kk < 0 || kk > nn) return -std::numeric_limits<double>::infinity();
    if (kk == 0 || kk == nn) return 0.0;
    return std::lgamma(nn + 1) - std::lgamma(kk + 1) - std::lgamma(nn - kk + 1);
}

// ranks[]: 1-based ranks (length n).  out[]: weights (length n), sum ~= 1.
inline void dnn_weights_from_ranks(const int *ranks, int n, int k, int s, double *out)
{
    validate_weight_parameters(n, k, s);
    validate_ranks(ranks, n);
    if (out == nullptr) throw std::invalid_argument("out must not be null");
    const double log_C_n_s = log_binom(n, s);
    double wsum = 0.0;

    for (int i = 0; i < n; ++i) {
        const int r = ranks[i];   // 1-based
        out[i] = 0.0;

        // Support: r <= n - s + k
        if (r > n - s + k) continue;

        // l bounds: l in [max(1, s-(n-r)), min(k, r, s)]
        int l_lo = std::max(1, s - (n - r));
        int l_hi = std::min({k, r, s});
        if (l_lo > l_hi) continue;

        double total = 0.0;
        for (int l = l_lo; l <= l_hi; ++l) {
            double lp = log_binom(r-1, l-1) + log_binom(n-r, s-l) - log_C_n_s;
            total += std::exp(lp);
        }
        out[i] = total / static_cast<double>(k);
        wsum  += out[i];
    }

    // Normalize for floating-point safety (should sum to 1 analytically)
    if (wsum > 0.0)
        for (int i = 0; i < n; ++i) out[i] /= wsum;
}

// ── Asymptotic DNN weights (Beta order-statistic increments) ──────────────
//
// w_rank[j] = (1/k) * Σ_{l=1}^{k} [F_Beta(l,s-l+1)(j/n) - F_Beta(l,s-l+1)((j-1)/n)]
//
// Full support over 1..n; normalized to sum 1.
inline void dnn_weights_asymp_from_ranks(const int *ranks, int n, int k, int s,
                                          double *out)
{
    validate_weight_parameters(n, k, s);
    validate_ranks(ranks, n);
    if (out == nullptr) throw std::invalid_argument("out must not be null");
    std::vector<double> w_rank(n, 0.0);

    for (int l = 1; l <= k; ++l) {
        double a = static_cast<double>(l);
        double b = static_cast<double>(s - l + 1);
        for (int j = 1; j <= n; ++j) {
            // F_Beta(a,b)(j/n) - F_Beta(a,b)((j-1)/n)
            // Use std::beta_distribution... not in standard C++. Use lgamma approach.
            // regularized incomplete beta: use Boost or manual.
            // For portability, compute via the log-beta PDF approximation.
            // ibeta(a, b, x) = regularized incomplete beta = I_x(a, b)
            // We use the digamma/lgamma series for the CDF difference.
            // Simple: CDF(j/n) - CDF((j-1)/n) ≈ pdf(j/n) * (1/n) for large n,
            // but we want it exact. Use recursion on the incomplete beta function.
            // Implementation: compute via the continued fraction or series.
            // For our use case (k,s moderate), use direct series evaluation.
            double x1 = static_cast<double>(j)     / n;
            double x0 = static_cast<double>(j - 1) / n;
            // regularized_ibeta via simple series (good for moderate a,b)
            auto reg_ibeta = [](double a_, double b_, double x) -> double {
                if (x <= 0.0) return 0.0;
                if (x >= 1.0) return 1.0;
                // Use continued fraction (Lentz) for stability
                // log B(a,b) = lgamma(a)+lgamma(b)-lgamma(a+b)
                double lbeta = std::lgamma(a_) + std::lgamma(b_) - std::lgamma(a_ + b_);
                // Use symmetry for x > (a_)/(a_+b_)
                bool sym = (x > (a_ / (a_ + b_)));
                if (sym) { std::swap(a_, b_); x = 1.0 - x; }
                // Lentz continued fraction for I_x(a,b)
                double qab = a_ + b_, qap = a_ + 1.0, qam = a_ - 1.0;
                double c = 1.0, d = 1.0 - qab * x / qap;
                if (std::fabs(d) < 1e-30) d = 1e-30;
                d = 1.0 / d;
                double h = d;
                for (int m = 1; m <= 200; ++m) {
                    double m2 = 2.0 * m;
                    double aa = m * (b_ - m) * x / ((qam + m2) * (a_ + m2));
                    d = 1.0 + aa * d; if (std::fabs(d) < 1e-30) d = 1e-30;
                    c = 1.0 + aa / c; if (std::fabs(c) < 1e-30) c = 1e-30;
                    d = 1.0 / d; h *= d * c;
                    aa = -(a_ + m) * (qab + m) * x / ((a_ + m2) * (qap + m2));
                    d = 1.0 + aa * d; if (std::fabs(d) < 1e-30) d = 1e-30;
                    c = 1.0 + aa / c; if (std::fabs(c) < 1e-30) c = 1e-30;
                    d = 1.0 / d;
                    double del = d * c;
                    h *= del;
                    if (std::fabs(del - 1.0) < 1e-12) break;
                }
                double result = std::exp(a_ * std::log(x) + b_ * std::log(1.0 - x)
                                         - lbeta) / a_ * h;
                return sym ? 1.0 - result : result;
            };
            w_rank[j - 1] += reg_ibeta(a, b, x1) - reg_ibeta(a, b, x0);
        }
    }

    // Divide by k, map to obs order, normalize
    double wsum = 0.0;
    for (int i = 0; i < n; ++i) {
        out[i] = w_rank[ranks[i] - 1] / k;
        wsum  += out[i];
    }
    if (wsum > 0.0)
        for (int i = 0; i < n; ++i) out[i] /= wsum;
}

// Full pipeline: X (n x p, column-major) + x0 (length p) -> weights (length n).
// inv_s: per-dimension inverse scales; nullptr for standard Euclidean.
// method: "exact" (default) or "asymp"
inline void dnn_weights_l2(const double *X, int n, int p,
                            const double *x0,
                            const double *inv_s,
                            int k, int s,
                            double *weights_out,
                            bool asymp = false)
{
    validate_sample_dimensions(n, p);
    validate_weight_parameters(n, k, s);
    validate_finite_array(X, static_cast<std::size_t>(n) * p, "X");
    validate_finite_array(x0, static_cast<std::size_t>(p), "x0");
    if (inv_s != nullptr) validate_finite_array(inv_s, static_cast<std::size_t>(p), "inv_s");
    if (weights_out == nullptr) throw std::invalid_argument("weights_out must not be null");
    std::vector<double> dist(n);
    std::vector<int>    ranks(n);
    dist_l2_core(X, n, p, x0, inv_s, dist.data(), /*squared=*/false);
    ranks_from_values(dist.data(), n, ranks.data());
    if (asymp)
        dnn_weights_asymp_from_ranks(ranks.data(), n, k, s, weights_out);
    else
        dnn_weights_from_ranks(ranks.data(), n, k, s, weights_out);
}

// Generic column-major n x (m*q) deletion matrix, in target/component blocks.
// Validation precedes output writes. Invert each permutation in O(n), then
// share the n-1 weight table across two O(n) passes per target/component.
inline void rank_loo(const int *ranks, const double *signal, const double *table,
                     int n, int m, int q, double *out)
{
    if (n<2 || m<1 || q<1) throw std::invalid_argument("rank_loo requires n>=2, m>=1 and q>=1");
    if (out==nullptr) throw std::invalid_argument("out must not be null");
    if (ranks==nullptr) throw std::invalid_argument("ranks must not be null");
    for (int j=0; j<m; ++j) validate_ranks(ranks+static_cast<std::size_t>(j)*n,n);
    validate_finite_array(signal,static_cast<std::size_t>(n)*q,"signal");
    validate_finite_array(table,static_cast<std::size_t>(n-1),"table");
    std::vector<int> ord(n);
    std::vector<double> left(n), right(n);
    for (int j=0; j<m; ++j) {
        const int *r=ranks+static_cast<std::size_t>(j)*n;
        for (int i=0; i<n; ++i) ord[r[i]-1]=i;
        for (int c=0; c<q; ++c) {
            const double *z=signal+static_cast<std::size_t>(c)*n;
            left[0]=0;
            for (int a=1; a<n; ++a) left[a]=left[a-1]+table[a-1]*z[ord[a-1]];
            right[n-1]=0;
            for (int a=n-2; a>=0; --a) right[a]=right[a+1]+table[a]*z[ord[a+1]];
            double *dest=out+(static_cast<std::size_t>(j)*q+c)*n;
            for (int i=0; i<n; ++i) dest[i]=left[r[i]-1]+right[r[i]-1];
        }
    }
}

} // namespace rdistnn
