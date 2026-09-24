# kernel/relaxation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/relaxation.m`
- Signature: `R=relaxation(spin_system,euler_angles)`
- Total lines: 897

## Purpose

Relaxation superoperator. Syntax: R=relaxation(spin_system,euler_angles)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

- SRSK contributes zero-destination relaxation rates to the accumulated generator; IME or DiBari-Levitt thermalisation is applied once to the accumulated spin relaxation, not separately to the recursive SRSK contribution. Bosonic mode dissipation is excluded from that recursive contribution and appended once after outer spin thermalisation, with the original mode parameters and temperature.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- euler_angles -three Euler angles (ZYZ active convention
- in radians) specifying system orientation
- relative to the input orientation; requi-
- by those theories that support relaxation
- rate anisotropy. It has no effect on tho-
- se theories (e.g. Redfield) that do not.

## Outputs

- R -relaxation superoperator. If a Liouvillian is
- assembled manually, this dissipative superoperator
- must enter as 1i*R, for example
- L=H+1i*R+1i*K; do not use H+R+K.
- Note: a variety of relaxation theories are supported, see the relax-
- ation theory parameters section of the online manual.
- Note: Spinach context functions include relaxation and kinetics
- superoperators into the total Liovillian automatically.
- Note: dissipative bosonic modes declared in inter.modes receive
- thermalised GKSL dissipators from rlx_modes.m in Liouville
- space formalisms; the euler_angles parameter refers to the
- spin subsystem only and has no effect on the mode terms.

## Header notes

Euler angles use the active ZYZ convention in radians and affect theories with anisotropic rates, not Redfield rotational averaging. In a manually assembled generator use L=H+1i*R+1i*K, not H+R+K.
