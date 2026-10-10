# tests/kernel/test_fourier_poly.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_fourier_poly.m)

`test_fourier_poly()` checks the production `fourdif` actions, adjoints, and transposes against explicit `fourdif` matrices for complex state blocks, odd/even grid sizes, derivative orders one to four, and non-unit periods. It then compares complete shaped RF pulse states and every trajectory point with explicit propagation across signed frequency slices, nonzero phase, multiple states, and two phase ranks. The `expv`, `evolution`, and `expm` methods are exercised, including the effective propagator retained by `expm` when polyadics are enabled. It returns the standard Spinach regression-result structure.
