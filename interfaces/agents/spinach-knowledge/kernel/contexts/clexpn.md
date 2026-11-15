# clexpn

Syntax: `[answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)`

The native context takes create output, an explicit registered five-argument experiment handle, native parameters and an interaction assumption string. Physical observed spin IDs live in parameters.observed; algorithm controls live in parameters.cce. No legacy input or result wrapper is accepted.

Closed-system stage: standard cluster expansion, mixed or explicit product bath states, Hilbert or Liouville generators and explicitly tagged projected branches. Dissipative, noise, sampled/interlaced context policies, partition expansions, native checkpoint/resume and density reduction are reserved for the next stage and fail explicitly here.
The context compiles immutable physical data, constructs retained generators, invokes the registered experiment and reduces raw responses. Natural numeric output is separate from diagnostics. Missing accuracy defaults are sqrt(eps); default memory budget is half the currently free physical memory. Missing order selects the included-bath order.
