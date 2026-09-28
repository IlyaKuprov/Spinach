# kernel/assume.m

- Signature: `spin_system=assume(spin_system,assumptions,retention)`

## Purpose

Selects the interaction-term approximations for a simulation context and stores them in `spin_system.inter.assumptions`. Call it before requesting the Hamiltonian.

## Parameters / inputs

- `spin_system` - spin-system structure to configure.
- `assumptions` - character string selecting the approximation set:
  - `nmr` - high-field NMR.
  - `esr` - electron rotating-frame ESR; `deer` uses the same set for DEER.
  - `deer-zz` - DEER with electron flip-flop terms removed.
  - `labframe` - full laboratory-frame Hamiltonian; bosonic modes and their interaction terms also remain in the laboratory frame.
  - `qnmr` - quadrupolar NMR with numerical rotating frames: spin-1/2 particles are in the rotating frame, while higher-spin particles start in the laboratory frame.
  - `cavity` - cavity QED: spins and bosonic modes share a rotating frame and the rotating-wave approximation. Mode energies are detunings from the carrier; exchange terms keep flip-flop components, while anharmonicity, Kerr, and dispersive terms are retained in full. Longitudinal and modulation terms are omitted as averaging out.
  - `spin-phonon` - spins use their usual rotating frames and bosonic modes stay in the laboratory frame. Electron and nuclear terms follow the ESR set; electron-mode exchange is omitted as non-secular, while mode-mode and nucleus-mode exchange, longitudinal, dispersive, modulation, and diagonal-mode terms are retained.
  - `se_dnp_h+`, `se_dnp_h-`, and `se_dnp_h0` - respectively select the positive-, negative-, and zero-frequency components of the solid-effect DNP Hamiltonian. They retain the corresponding `EzNp`/`EzNm`/`EzNz` electron-nuclear terms and `T(L,+1)`/`T(L,-1)`/secular inter-nuclear terms; inter-electron, giant-spin, quadratic, and Zeeman interactions are ignored.
- `retention` - optional spin-only retention mask. `zeeman` drops spin-spin interactions; `couplings` drops Zeeman interactions. These retention options are undefined for systems containing bosonic modes.

## Outputs

- Updated `spin_system`, with the selected interaction strengths configured.

## Reference

[Spin Dynamics Wiki: assume.m](https://spindynamics.org/wiki/index.php?title=assume.m)
