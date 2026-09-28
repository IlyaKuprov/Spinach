# kernel/create.m

- Signature: `spin_system=create(sys,inter)`

## Purpose

The entry function of the Spinach kernel that constructs the spin-system object required by the rest of the library. It validates and absorbs the spin-system, instrument, and interaction specifications, then writes diagnostics to the console.

## Parameters / inputs

- `sys` — spin-system and instrument specification structure; see the spin-system specification section of the online manual.
- `inter` — interaction specification structure; see the spin-system specification section of the online manual.

## Output

- `spin_system` — the primary object used by Spinach to store simulation information.

## Bosonic mode notes

- `inter.modes.carriers` declares, for each bosonic mode, the laboratory frequency of the rotating frame in which `inter.modes.frqs` is specified. The declared frequencies are detunings and may be negative; thermal occupations are computed using the physical frequency, the sum of the carrier and the detuning.
- `inter.modes.t2_times` is interpreted at the declared temperature. The pure-dephasing rate is `1/T2-kappa*(1+2*nbar)/2`, where `kappa` is the amplitude damping rate and `nbar` is the thermal occupation at the physical mode frequency.
- Bosonic-mode quadrature operators use the normalization `(a+a')/sqrt(2)` in both `inter.modes.longitudinal` and the modulation channels.
