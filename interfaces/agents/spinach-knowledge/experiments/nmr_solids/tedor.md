# experiments/nmr_solids/tedor.m

Use with Liouville-space `singlerot` for a pair of distinct spin-1/2 isotopes. The direct sequence builds antiphase coherence on the first channel, converts with simultaneous finite quarter turns, and detects the second channel immediately and then once per rotor period using the sender inversion train.

Both zero-filter (ZF) sequences detect the prepared channel after two recoupling sandwiches. `zf-cycle` subtracts the opposite-phase branch from the equal-phase branch with the receiver signs included. `zf-filter` instead retains the first-channel single-quantum coherences that are diagonal on the partner. They are not interchangeable. Complex outputs are kept through powder averaging.
