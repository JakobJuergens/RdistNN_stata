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
//
// Ported from RdistNN/src/{dist_l2.cpp, fast_order.cpp} with Rcpp removed.

#pragma once

#include <vector>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <limits>

namespace rdistnn {

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
// Support: r = 1, …, n-s+k  (NOT 1..s).
// For k=1: w_r = C(n-r, s-1) / C(n, s)  (standard Steele formula).
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
    std::vector<double> dist(n);
    std::vector<int>    ranks(n);
    dist_l2_core(X, n, p, x0, inv_s, dist.data(), /*squared=*/false);
    ranks_from_values(dist.data(), n, ranks.data());
    if (asymp)
        dnn_weights_asymp_from_ranks(ranks.data(), n, k, s, weights_out);
    else
        dnn_weights_from_ranks(ranks.data(), n, k, s, weights_out);
}

} // namespace rdistnn
