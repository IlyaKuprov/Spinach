# kernel/homospoil.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/homospoil.m`
- Signature: `rho=homospoil(spin_system,rho,zqc_flag)`
- Total lines: 118

## Purpose

Emulates a strong homospoil pulse -only zero-frequency states with respect to the carrier frequencies (chemical shifts are not conside- red) survive the process. Syntax: rho=homospoil(spin_system,rho,zqc_flag)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- rho -a state vector or a horizontal stack thereof
- zqc_flag -a flag controlling the fate of zero-quantum
- coherences. If set to 'keep', causes ZQCs to
- survive the process, approximating experimen-
- tal behaviour. If set to 'destroy', wipes the
- zero-quantum coherences -only the longitudi-
- nal states survive the process.
- The flag is ignored in zeeman-hilb and zeeman-
- liouv formalisms, where the effect is always
- to destroy everything except the diagonal of
- the density matrix.

## Outputs

- rho -the state vector(s) with only the longitudi-
- nal or only the zero-quantum states kept
- Note: this function is only available for sphten-liouv formalism; it
- supports Fokker-Planck direct products.
- Note: this is a purely mathematical filter that only mimics -in an
- idealised way -the effect of a real homospoil pulse. Essenti-
- ally, it searches the density matrix for any transverse state
- populations and zeroes them out. If the flag is set, zero-qua-
- ntum coherences are also erased.

## Implementation structure

- Emulates a strong homospoil pulse -only zero-frequency states with
- respect to the carrier frequencies (chemical shifts are not conside-
- red) survive the process. Syntax:
- rho=homospoil(spin_system,rho,zqc_flag)
- rho -a state vector or a horizontal stack thereof
- zqc_flag -a flag controlling the fate of zero-quantum
- coherences. If set to 'keep', causes ZQCs to
- survive the process, approximating experimen-
- tal behaviour. If set to 'destroy', wipes the
- zero-quantum coherences -only the longitudi-
- nal states survive the process.
- The flag is ignored in zeeman-hilb and zeeman-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `spdiags()`, `true()`, `offdiag_mask()`, `rho()`, `lin2lm()`, `report()`, `vector()`, `ischar()`, `ismember()`.
