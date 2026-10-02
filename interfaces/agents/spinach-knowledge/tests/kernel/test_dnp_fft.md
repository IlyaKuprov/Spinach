# tests/kernel/test_dnp_fft.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dnp_fft.m)

`test_dnp_fft()` checks: Complete electron-proton steady-state frequency scans in both Liouville bases, signed offsets, multiple coils, and odd/even phase grids. Implicit GMRES is compared with explicit backslash and the original matrix GMRES; direct FP and LvN behaviour is checked with the polyadic enable switch on and off. It returns the standard Spinach regression-result structure.
