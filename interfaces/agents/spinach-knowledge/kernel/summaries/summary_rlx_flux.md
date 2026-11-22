# kernel/summaries/summary_rlx_flux.m

`summary_rlx_flux(spin_system)` reports validated point-to-point magnetisation flux rates, direction and correlation retention. It reports the source, destination, and rate from each compiled record, including diagonal records. `create` calls it directly before basis construction. Absent or zero flux prints nothing. No generator is constructed or changed.

Flux reporting is gated by selection of 'flux' in spin_system.rlx.theories; disabled flux is not reported merely because rate fields exist. Active rates and types live in spin_system.rlx.flux_rate and spin_system.rlx.flux_type, not chemistry. This routine reports only and never constructs or appends the operator.
