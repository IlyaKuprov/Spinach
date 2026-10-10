# tests/kernel/test_sim2liouv_cache.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_sim2liouv_cache.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_sim2liouv_cache.m)

## Purpose

Regression test for the Hilbert-to-Liouville conversion function `sim2liouv`, focused on cache identity across representations. The test verifies that converted operators and Hamiltonians use the converted basis cache identity, that cache insertion order does not cause collisions, that cache metadata remains correct when caching is disabled, that objects which never requested a cache are preserved, and that cache use does not change acquired signals.

## Behaviour

- Registers a test named `kernel/sim2liouv_cache` with the description "Converted basis cache identity".
- Builds small single-spin systems (`sys.magnet=1`, `sys.isotopes={'1H'}` or `{'13C'}`, `sys.output='hush'`, `sys.disable={'hygiene'}`, `sys.parallel={'processes',1}`, `sys.parprops={}`) with `inter.zeeman.scalar={1}` and `bas.formalism='zeeman-hilb'`, `bas.approximation={'none'}`.
- Iterates over cache flag combinations `{'op_cache'}`, `{'ham_cache'}`, and `{'op_cache','ham_cache'}`, and over two cache-warming orderings (Hilbert first or Liouville first), setting `inter.zeeman.scalar={flag_idx}` and `sys.isotopes=nuclei(ordering)` to separate fixtures without altering physical cache keys.
- For each combination, converts the Hilbert-space system to Liouville space via `sim2liouv(spin_h,struct(),H,[],[])`, and checks that the converted basis hash equals `md5_hash({spin_l.bas.basis,spin_l.bas.nstates,spin_l.chem.parts})` ("canonical converted hash").
- Warms both representations in the requested order using `operator(...,'L+',1)` and `hamiltonian(...)`, then compares the cached raising operator against the reference operator and the cached drift Hamiltonian against the reference Hamiltonian with tolerances `1e-12` (absolute and relative). The reference for ordering 1 is the Liouville operator and the conversion map `K`; for ordering 2 it is the Hilbert operator and Hamiltonian `H`.
- With caching disabled (`sys.enable={}`), converts again with empty Hamiltonian inputs and checks that the existing basis hash still matches `md5_hash({spin_l.bas.basis,spin_l.bas.nstates,spin_l.chem.parts})` ("disabled cache hash"), so metadata remains valid if caching is re-enabled.
- For a system that never requested a cache, checks the same canonical compiled hash ("uncached canonical hash"); identity is independent of whether caching is enabled.
- For the three no-op formalisms `zeeman-liouv`, `sphten-liouv`, and `zeeman-wavef`, calls `sim2liouv(spin_system,parameters,H,R,K)` with `H=operator(spin_system,'Lx',1)+0.37*operator(spin_system,'Ly',1)`, `R` and `K` empty sparse matrices, and `parameters.tag='unchanged'`, verifying that the spin system, parameters, Hamiltonian, relaxation matrix, and auxiliary matrix are returned unchanged ("no-op inputs").
- Acquires a complex FID after a noncommuting phase-shifted soft pulse using `liquid(spin_system,@sp_acquire,parameters,'nmr')` with `sys.enable={'op_cache','ham_cache'}`, `inter.zeeman.scalar={0}`, `parameters.spins={'1H'}`, `parameters.rho0=state(spin_system,'Lz','1H')`, `parameters.coil=state(spin_system,'L+','1H')`, `parameters.pulse_frq=17`, `parameters.pulse_phi=0.37`, `parameters.pulse_pwr=2*pi*100`, `parameters.pulse_dur=0.0025`, `parameters.pulse_rnk=2`, `parameters.offset=31`, `parameters.sweep=1000`, `parameters.npoints=8`, `parameters.method='expm'`. Compares the cached-enable FID against the cache-disabled FID with tolerances `1e-11` ("cached soft-pulse FID"), and checks `norm(fid_ref)>0.1` and `norm(imag(fid_ref))>0.1` ("nonzero complex signal").

## Inputs and outputs

```matlab
result = test_sim2liouv_cache()
```

- **Inputs:** none.
- **Outputs:** `result` — a test result object accumulated via `test_close`, containing regression checks for both cache insertion orders, complex operators, soft-pulse acquisition, and no-op paths.

## References

- [sim2liouv](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sim2liouv.m) — Hilbert-to-Liouville conversion under test.
- [Spinach GitHub repository](https://github.com/IlyaKuprov/Spinach)

A segmented two-plus-four-dimensional fixture rejects nonzero entries in both off-diagonal substance blocks of H, R, K, all four operator-like fields, and the second member of each horizontal state stack. A valid sparse block-diagonal generator, pulse operator, and two-state stack are compared exactly with independent substance conversions.
