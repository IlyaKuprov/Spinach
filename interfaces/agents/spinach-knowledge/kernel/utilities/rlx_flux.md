# kernel/utilities/rlx_flux.m

`F=rlx_flux(spin_system)` constructs the point-to-point magnetisation flux relaxation superoperator from validated `chem.flux_rate` and `chem.flux_type` fields and compiled local basis descriptors. Rows of the rate matrix are sources and columns are destinations.

The zero-flux path returns an allocated zero generator. Nonzero flux requires `sphten-liouv`. Intramolecular flux preserves transferred correlations; intermolecular flux drains source correlations and transfers single-spin terms. Direct-sum blocks use local spin membership and offsets. `relaxation` adds this contribution after spin thermalisation and bosonic dissipation, or returns it alone when ordinary relaxation is disabled.
