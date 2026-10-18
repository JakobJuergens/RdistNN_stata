#include "rdistnn_core.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
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

template <typename Fn>
bool rejects(Fn fn)
{
    try {
        fn();
    } catch (const std::invalid_argument &) {
        return true;
    }
    return false;
}

} // namespace

int main()
{
    // Two targets and paired signals: catch shifted-suffix and scatter errors.
    const int jk_ranks[] = {1,2,3,4, 4,3,2,1};
    const double jk_signals[] = {2,1,4,3, 1,-1,2,1};
    const double jk_weights[] = {2.0/3,1.0/3,0};
    const double jk_expected[] = {
        2,8.0/3,5.0/3,5.0/3, 0,4.0/3,1.0/3,1.0/3,
        10.0/3,10.0/3,7.0/3,3, 4.0/3,4.0/3,1.0/3,1
    };
    std::vector<double> jk_out(16);
    rdistnn::rank_loo(jk_ranks, jk_signals, jk_weights, 4, 2, 2, jk_out.data());
    for (std::size_t i=0; i<jk_out.size(); ++i)
        if (!close(jk_out[i], jk_expected[i])) return fail("prefix/suffix paired deletions");
    const int bad_jk_ranks[] = {1,1,3,4};
    std::fill(jk_out.begin(),jk_out.end(),17.0);
    if (!rejects([&] { rdistnn::rank_loo(bad_jk_ranks,jk_signals,jk_weights,4,1,2,jk_out.data()); }))
        return fail("jackknife accepted duplicate ranks");
    for (double value: jk_out)
        if (value!=17.0) return fail("invalid jackknife modified output");

    std::vector<int> support_ranks(30), support_ranks_m1(29);
    std::iota(support_ranks.begin(), support_ranks.end(), 1);
    std::iota(support_ranks_m1.begin(), support_ranks_m1.end(), 1);
    std::vector<double> support_signal(30, 1e16), support_table(29), support_out(30);
    support_signal.back() = -1e307;
    rdistnn::dnn_weights_from_ranks(support_ranks_m1.data(),29,1,10,support_table.data());
    rdistnn::rank_loo(support_ranks.data(),support_signal.data(),support_table.data(),
                     30,1,1,support_out.data());
    for (double value: support_out)
        if (value!=support_out[0]) return fail("unused outlier creates native deletion variation");
    // Subtracting the reference from the unused last row would overflow.
    std::fill(support_signal.begin(),support_signal.end(),1e308);
    support_signal.back()=-1e308;
    rdistnn::rank_loo(support_ranks.data(),support_signal.data(),support_table.data(),
                     30,1,1,support_out.data());
    for (double value: support_out)
        if (!std::isfinite(value) || value!=support_out[0])
            return fail("unused extreme overflows native centering");

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

    const std::vector<double> sentinel(n, 17.0);
    exact = sentinel;
    if (!rejects([&] { rdistnn::dnn_weights_from_ranks(ranks.data(), n, 0, 2, exact.data()); }))
        return fail("k=0 was not rejected");
    if (exact != sentinel) return fail("k=0 modified output");
    if (!rejects([&] { rdistnn::dnn_weights_from_ranks(ranks.data(), n, 3, 2, exact.data()); }))
        return fail("k>s was not rejected");
    if (!rejects([&] { rdistnn::dnn_weights_from_ranks(ranks.data(), n, n + 1, n + 1, exact.data()); }))
        return fail("k>n was not rejected");
    if (!rejects([&] { rdistnn::dnn_weights_from_ranks(ranks.data(), n, 1, n + 1, exact.data()); }))
        return fail("s=n+1 was not rejected");
    if (!rejects([&] { rdistnn::dnn_weights_from_ranks(ranks.data(), 0, 1, 1, exact.data()); }))
        return fail("empty sample was not rejected");
    if (!rejects([&] { rdistnn::dnn_weights_l2(X, n, 0, x0, nullptr, 1, 1, exact.data(), false); }))
        return fail("p=0 was not rejected");
    const double bad_x0[] = {std::numeric_limits<double>::quiet_NaN()};
    if (!rejects([&] { rdistnn::dnn_weights_l2(X, n, p, bad_x0, nullptr, 1, 1, exact.data(), false); }))
        return fail("nonfinite x0 was not rejected");

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
