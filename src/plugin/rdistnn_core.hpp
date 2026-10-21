// Frozen-rank deletion kernel for the RdistNN Stata plugin.
// Distances live in rdistnn_distance.hpp; Mata supplies ranks and weight tables.
#pragma once

#include "compensated_arithmetic.hpp"
#include "rdistnn_distance.hpp"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace rdistnn {

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
    rdistnn_arithmetic::rank_loo(ranks,signal,table,n,m,q,out);
}

} // namespace rdistnn
