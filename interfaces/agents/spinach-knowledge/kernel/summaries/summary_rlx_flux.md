# kernel/summaries/summary_rlx_flux.m

`summary_rlx_flux(spin_system)` reports validated point-to-point magnetisation flux rates, direction and correlation retention. Rows are sources and columns destinations. `create` calls it directly before basis construction. Absent or zero flux prints nothing. No generator is constructed or changed.
