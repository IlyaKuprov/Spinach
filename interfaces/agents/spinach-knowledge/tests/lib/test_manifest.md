# tests/lib/test_manifest.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

The regression registry maps stable kernel test IDs and descriptions to test functions. The existing CCE-family entries are retained, with kernel/clexpn_options and kernel/clexpn_transition added for the staged native input services. The new functions return standard result records and are executable through run_test and run_tests; they do not replace or weaken the existing CCE assertions.
