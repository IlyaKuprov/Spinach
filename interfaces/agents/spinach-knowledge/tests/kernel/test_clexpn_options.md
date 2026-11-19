# tests/kernel/test_clexpn_options.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_clexpn_options.m)

Registered as kernel/clexpn_options. Repeats the finite-bath lattice, occupancy, frame, total-Hz-tensor, density, isotopologue, exclusion, group, and seeded-stream assertions through spin_bath. The new input contract additionally checks declared defaults, explicit populations and jump rates, rejected old ownership and unknown fields, invalid densities/IDs/seeds, explicit pump frequency units and species filters, unchanged global random state, preservation of native properties on append, and generalised physical densities without an unnecessary transition selector. The column-isotope fixture is explicitly supplied as the documented cell row; production code does not normalise input shape.

Native `bath_jumps` retains the complete physical GKSL support contract: selected observed spins as well as selected bath spins are allowed. Generalised dynamics may transfer observed populations; projected dynamics still rejects jumps that mix the declared manifolds. IDs outside the selected observed-plus-bath set are rejected.
