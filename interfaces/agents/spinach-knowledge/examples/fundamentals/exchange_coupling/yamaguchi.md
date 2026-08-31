# examples/fundamentals/exchange_coupling/yamaguchi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/exchange_coupling/yamaguchi.m`
- Signature: `yamaguchi()`
- Total lines: 18

## Purpose

Yamaguchi equation estimate of exchange coupling from a broken-symmetry DFT calculation on a bistrityl bira- dical with an alkynyl linker from Olav Schiemann.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Read Gaussian logs; implemented by `props_sing=gparse('biradical_singlet.log')`.
- Lines 13-14: Call Yamaguchi equation; implemented by `J=brokensymm(props_sing,props_trip)`.

### Key state/data transformations

- Lines 10: computes `props_sing` using `props_sing=gparse('biradical_singlet.log')`.
- Lines 11: computes `props_trip` using `props_trip=gparse('biradical_triplet.log')`.
- Lines 14: computes `J` using `J=brokensymm(props_sing,props_trip)`.

## Implementation structure

- Yamaguchi equation estimate of exchange coupling from
- a broken-symmetry DFT calculation on a bistrityl bira-
- dical with an alkynyl linker from Olav Schiemann.
- Read Gaussian logs
- Call Yamaguchi equation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `brokensymm()`, `num2str()`.
