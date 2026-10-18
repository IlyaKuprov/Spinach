# kernel/summaries/summary_kinetics.m

`summary_kinetics(spin_system)` reports declared chemistry, matched spins, traced spins and closure directly from validated `chem.reactions` and `chem.parts`. `create` calls it before basis construction; no compiled reaction operators are needed. Flux reporting remains separate in `summary_rlx_flux`.
