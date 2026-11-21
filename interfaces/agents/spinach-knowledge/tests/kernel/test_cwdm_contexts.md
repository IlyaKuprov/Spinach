# tests/kernel/test_cwdm_contexts.m

Calls all nine contexts that consume `kinetics` matrices (`liquid`, `imaging`, `crystal`, `powder`, `device`, `floquet`, `singlerot`, `doublerot`, and `gridfree`) with callback-rate and bimolecular reaction records. Checks the context-specific error identifier and the actionable custom-sequence guidance. Constant first-order kinetics must still reach the supplied sequence. Small rotor and spatial grids test assembly, not acquisition accuracy.

The imaging cases select `'flux'` in `inter.relaxation` and use top-level `inter.flux_rate` and `inter.flux_type` to construct a relaxation operator. With empty phantom operators, imaging returns zero relaxation despite active flux metadata. When a complete operator containing flux is supplied, imaging applies exactly the caller's spatial weights and adds no flux to either `R` or `K`. Zero-rate Zeeman input remains valid.
