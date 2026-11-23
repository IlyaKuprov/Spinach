# periodic_acquire.m

The Hilbert-space single-rotor callback evaluates coherent periodic FIDs using a one-period unitary and within-period prefixes. Dwell must be an integer number of ticks. Rotor Hamiltonians are held from their Fourier node to the next tick, and start-tick changes cyclically shift RF and drift together. Nonzero relaxation and kinetics are rejected.
