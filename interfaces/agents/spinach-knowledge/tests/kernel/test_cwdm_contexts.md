# tests/kernel/test_cwdm_contexts.m

Calls all nine contexts that consume `kinetics` matrices (`liquid`, `imaging`, `crystal`, `powder`, `device`, `floquet`, `singlerot`, `doublerot`, and `gridfree`) with callback-rate and bimolecular reaction records. Checks the context-specific error identifier and the actionable custom-sequence guidance. Constant first-order kinetics must still reach the supplied sequence. Small rotor and spatial grids test assembly, not acquisition accuracy.
