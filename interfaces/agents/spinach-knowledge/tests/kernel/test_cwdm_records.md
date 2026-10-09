# tests/kernel/test_cwdm_records.m

Validates explicit reaction input, the additive default, retirement errors, matching ownership/isotopes/uniqueness, scalar or time-dependent rates, selector structure, and merge offsets. Every new create rejection is asserted by identifier and message; accepted loss, spin-free sink, repeated-substance, and selector records exercise the production parser and summary. This is an input-contract test, not a generator-physics test.
