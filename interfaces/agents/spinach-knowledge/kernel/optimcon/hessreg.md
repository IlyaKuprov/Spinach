# kernel/optimcon/hessreg.m

- Signature: `[H,data]=hessreg(spin_system,H,g,data)`

## Purpose

RFO regularisation for Newton-Raphson Hessian and gradient pairs. Syntax: [H,data]=hessreg(spin_system,H,g,data)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

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
