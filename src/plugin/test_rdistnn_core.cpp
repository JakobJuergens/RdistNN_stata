#include "rdistnn_core.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

bool close(double actual, double expected, double tolerance = 1e-12)
{
    return std::fabs(actual - expected) <= tolerance;
}

int fail(const char *message)
{
    std::cerr << "rdistnn core test failed: " << message << '\n';
    return 1;
}

} // namespace

int main()
{
    const int n = 4;
    const int p = 1;
    const double X[] = {0.0, 1.0, 2.0, 3.0};
    const double x0[] = {0.0};
    std::vector<double> distances(n);
    std::vector<int> ranks(n);

    rdistnn::dist_l2_core(X, n, p, x0, nullptr, distances.data());
    for (int i = 0; i < n; ++i) {
        if (!close(distances[i], static_cast<double>(i)))
            return fail("L2 distance mismatch");
    }

    rdistnn::ranks_from_values(distances.data(), n, ranks.data());
    for (int i = 0; i < n; ++i) {
        if (ranks[i] != i + 1)
            return fail("rank mismatch");
    }

    std::vector<double> exact(n);
    rdistnn::dnn_weights_from_ranks(ranks.data(), n, 1, n, exact.data());
    if (!close(exact[0], 1.0))
        return fail("s=n endpoint did not select the nearest observation");
    for (int i = 1; i < n; ++i) {
        if (!close(exact[i], 0.0))
            return fail("s=n endpoint has weight outside the nearest rank");
    }

    rdistnn::dnn_weights_from_ranks(ranks.data(), n, 2, 2, exact.data());
    for (double weight : exact) {
        if (!close(weight, 0.25))
            return fail("k=s did not produce uniform weights");
    }

    std::vector<double> asymp(n);
    rdistnn::dnn_weights_asymp_from_ranks(ranks.data(), n, 1, 2, asymp.data());
    const double asymp_sum = std::accumulate(asymp.begin(), asymp.end(), 0.0);
    if (!close(asymp_sum, 1.0))
        return fail("asymptotic weights are not normalized");
    for (double weight : asymp) {
        if (!std::isfinite(weight) || weight < 0.0)
            return fail("asymptotic weights are invalid");
    }

    const double tied[] = {1.0, 1.0, 0.0, 2.0};
    const int expected_tied_ranks[] = {2, 3, 1, 4};
    rdistnn::ranks_from_values(tied, n, ranks.data());
    for (int i = 0; i < n; ++i) {
        if (ranks[i] != expected_tied_ranks[i])
            return fail("stable tie ordering mismatch");
    }

    std::vector<double> pipeline(n);
    rdistnn::dnn_weights_l2(
        X, n, p, x0, nullptr, 1, n, pipeline.data(), false
    );
    for (int i = 0; i < n; ++i) {
        const double expected = (i == 0) ? 1.0 : 0.0;
        if (!close(pipeline[i], expected))
            return fail("full L2 pipeline mismatch");
    }

    std::cout << "rdistnn core test passed\n";
    return 0;
}
