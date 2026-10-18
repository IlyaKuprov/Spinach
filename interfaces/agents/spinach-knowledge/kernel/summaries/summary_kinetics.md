# kernel/summaries/summary_kinetics.m

`summary_kinetics(spin_system,reactions)` reports the declared chemistry and each compiled reaction’s matched spins, traced spins, and closure. Both arguments are required. `kinetics(spin_system,'report')` calls it after reaction compilation. Reporting does not construct or alter the generator. Flux reporting is separate in `summary_rlx_flux`.
