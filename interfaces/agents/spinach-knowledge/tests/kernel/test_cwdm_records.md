# tests/kernel/test_cwdm_records.m

Validates explicit reaction input, the additive default, retirement errors, matching ownership/isotopes/uniqueness, scalar or time-dependent rates, selector structure, and merge offsets. Every new create rejection is asserted by identifier and message; accepted loss, spin-free sink, repeated-substance, and selector records exercise the production parser and summary. This is an input-contract test, not a generator-physics test.

An empty-chemistry merge regression checks that absent partitions and reaction records remain absent and that `create` supplies one default substance containing the merged spins.

Column-oriented single-substance membership is preserved by create and printed on one line by the chemistry summary.
