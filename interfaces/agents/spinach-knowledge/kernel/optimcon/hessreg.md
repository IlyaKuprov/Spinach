# kernel/optimcon/hessreg.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/hessreg.m`
- Signature: `[H,data]=hessreg(spin_system,H,g,data)`
- Total lines: 101

## Purpose

RFO regularisation for Newton-Raphson Hessian and gradient pairs. Syntax: [H,data]=hessreg(spin_system,H,g,data)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- H -Hessian matrix to be regularised
- g -gradient computed at the same point as H
- data -diagnostic data structure

## Outputs

- H -regularised Hessian
- data -updated diagnostic data structure with
- data.count.rfo incremented by the number
- of RFO iterations taken

## Implementation structure

- RFO regularisation for Newton-Raphson Hessian and gradient
- pairs. Syntax:
- [H,data]=hessreg(spin_system,H,g,data)
- H -Hessian matrix to be regularised
- g -gradient computed at the same point as H
- data -diagnostic data structure
- H -regularised Hessian
- data -updated diagnostic data structure with
- data.count.rfo incremented by the number
- of RFO iterations taken
- Check consistency
- Set shorthands

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `logical()`, `cond()`, `speye()`, `issymmetric()`, `iscolumn()`.
