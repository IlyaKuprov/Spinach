# kernel/summaries/summary_chemistry.m

`summary_chemistry(spin_system)` reports each substance's global spin membership and initial concentration, including spin-free substances. Every explicit reaction is printed with reactant and product substance indices, scalar rate or time-handle expression, closure, matching pairs, and optional selector name (or user product-superoperator pair). Column-oriented spin memberships are formatted as rows without mutating their stored shape. Output is routed through `report`; the system is not modified. Empty products denote loss without a tracked product.
