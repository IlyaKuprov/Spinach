# kernel/summaries/summary_kinetics.m

`summary_kinetics(spin_system)` reports validated point-to-point magnetisation flux rates, their direction, and whether correlations are retained. It works before basis construction and prints nothing for absent or zero flux.

`summary_kinetics(spin_system,reactions)` reports the declared chemistry and each compiled reaction’s matched spins, traced spins, and closure. `kinetics(spin_system,'report')` calls these forms before and after reaction compilation, respectively. Reporting does not construct or alter the generator.
