// Public CATE/NP plot order must not depend on which ado was autoloaded first.
version 17
do "tests/test_framework.do"

test_section "plots: [autoload] CATE before NP in a fresh runner session"
rdnn_dgp_basic, n(80) seed(52)
quietly cate_dnn y, xcols(x1 x2) x0(0 0) ///
    tcol(D) ecol(e) m0col(m0) m1col(m1)
capture noisily cate_dnn_plot, type(estimates) name(order_cate_first) replace nodraw
local cate_first_rc = _rc
test_assert "[autoload] first CATE plot builds" `cate_first_rc' == 0

quietly np_dnn y, xcols(x1 x2) x0(0 0)
matrix order_np_b = e(b)
capture noisily np_dnn_plot, type(estimates) name(order_np_next) replace nodraw
local np_next_rc = _rc
local np_next_n = r(n_rows)
test_assert "[autoload] NP plot builds after CATE" `np_next_rc' == 0
test_assert "[autoload] NP plot reports its target" `np_next_n' == 1
test_assert "[autoload] NP plot preserves estimates" _b[mu] == order_np_b[1,1]

quietly cate_dnn y, xcols(x1 x2) x0(0 0) ///
    tcol(D) ecol(e) m0col(m0) m1col(m1)
capture noisily cate_dnn_plot, type(estimates) name(order_cate_again) replace nodraw
local cate_again_rc = _rc
test_assert "[autoload] CATE plot still builds after NP" `cate_again_rc' == 0
