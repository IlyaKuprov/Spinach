# kernel/states/four_spin_states.m

- MATLAB implementation: [kernel/states/four_spin_states.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/four_spin_states.m)

## Signature

`rho = four_spin_states(spin_system, spins, spin_state)`

## Purpose

Constructs one of sixteen singlet/triplet product states for four selected spin-1/2 particles. The selector assigns a state to each of two ordered pairs. The function assembles the result from explicit Cartesian spin-operator expansions in the active Spinach representation; it constructs a state and does not propagate it. The MATLAB source header also points to an accompanying Mathematica file.

## Inputs

- `spin_system` — a Spinach system used by `state` and by the consistency checks on spin count and multiplicities.
- `spins` — a 1-by-4 row vector of four distinct positive real integer spin indices. Every index must be at most `spin_system.comp.nspins`, and each indexed spin must have multiplicity 2 (spin-1/2). These are indices, not spin quantum numbers. Their supplied order defines the pairs: the first pair is `spins(1:2)`; the second is `spins(3:4)`.
- `spin_state` — a MATLAB character value equal to one of the exact, case-sensitive selectors in the table. The source rejects non-character input and errors on an unrecognised selector.

## Supported pair-state selectors

In each selector, the token before `(x)` describes the first pair `spins(1:2)`; the token after it describes the second pair `spins(3:4)`. `S` means singlet, `TU` triplet-up, `T0` triplet-middle, and `TD` triplet-down.

| First pair \ Second pair | `S` | `TU` | `T0` | `TD` |
|---|---|---|---|---|
| `S`  | `S(x)S`  | `S(x)TU`  | `S(x)T0`  | `S(x)TD` |
| `TU` | `TU(x)S` | `TU(x)TU` | `TU(x)T0` | `TU(x)TD` |
| `T0` | `T0(x)S` | `T0(x)TU` | `T0(x)T0` | `T0(x)TD` |
| `TD` | `TD(x)S` | `TD(x)TU` | `TD(x)T0` | `TD(x)TD` |

The selector is ordered: reversing the two pair states requires using the correspondingly reversed selector and pair ordering.

## Construction and output

The routine calls `grumble` to enforce the input constraints, then builds one-, two-, three-, and four-spin operator components using `state(spin_system, ... )` with `Lx`, `Ly`, `Lz`, and identity `E` factors assigned in the order of `spins(:)`. A selector-specific explicit expansion of those components is returned as `rho`.

The source documents `rho` as a density matrix in Hilbert formalism or a state vector in Liouville formalism. Its representation is therefore the one produced through Spinach's `state` machinery for the supplied system; the kernel does not promise one fixed array shape independent of formalism and basis. The available choices are limited to the sixteen pair-product states listed above, and all four selected spins must be spin-1/2.
