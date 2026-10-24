#ifndef RDISTNN_TEST_STATA_SDK_HPP
#define RDISTNN_TEST_STATA_SDK_HPP

// Independent test-host model, calibrated against the licensed Stata SDK.
// Do not reuse a production range constant: both sides could then share a bug.
namespace rdistnn_test_sdk {
inline bool is_missing(double value) { return value >= 0x1p1023; }
}

#endif
