# tests/kernel/test_prop_unrounded.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_prop_unrounded.m)

`result=test_prop_unrounded()` checks sparse and dense nonnormal exponentials without chopping, using disabled cleanup or a zero chop tolerance, as well as signed/complex timesteps, a nilpotent control, the normal chopped path, and a shared-cache miss/hit. It uses the standard regression result helpers and an explicitly sized Spinach process pool. The unrounded results are compared with MATLAB expm; the test must be run under the normal test runner wall-clock budget because a termination regression can stall.
