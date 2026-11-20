# Applying a Cartesian dipolar symbol

Prepare the spatial multiplier once with `ddf_bounded` or `ddf_periodic`. At each dynamical stage call `ddf_cartesian(magn,npts,symbol)` with physical magnetisation in A/m and X varying fastest. The tensor is traceless after secular averaging: the transverse multipliers are each minus half the longitudinal one. All three output columns are in tesla; target-isotope gyromagnetic ratios belong in the spin dynamics, not in the field calculation.
