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

double choose_small(int n, int k)
{
    unsigned long long value=1;
    for (int j=1; j<=k; ++j) value=value*(n-k+j)/j;
    return static_cast<double>(value);
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

    std::vector<int> support_ranks(30);
    std::iota(support_ranks.begin(), support_ranks.end(), 1);
    std::vector<double> support_signal(30, 1e16), support_table(29), support_out(30);
    support_signal.back() = -1e307;
    // Independent small-integer subset probabilities; table construction is
    // owned by Mata, so this native deletion fixture supplies its own table.
    for (int i=0; i<20; ++i)
        support_table[i]=choose_small(28-i,9)/choose_small(29,10);
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

    // Independent binary-power anchor: deleting the zero leaves A+1-A=1.
    const int cancellation_ranks[] = {1, 2, 3, 4};
    const double A = std::ldexp(1.0, 54);
    const double cancellation_signal[] = {A, 1, 0, -A, 7, 7, 7, 7};
    const double cancellation_table[] = {1, 1, 1};
    double cancellation_loo[8];
    rdistnn::rank_loo(cancellation_ranks, cancellation_signal,
        cancellation_table, 4, 1, 2, cancellation_loo);
    if (cancellation_loo[2] != 1 || cancellation_loo[6] != 21)
        return fail("prefix/suffix merge discarded a correction");

    std::cout << "rdistnn core test passed\n";
    return 0;
}
