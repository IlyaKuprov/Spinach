# kernel/optimcon/wrappers/grape_phase.m

- Signature: `[traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)`

## Purpose

Cost function for optimal control using the GRAPE algorithm. Returns fidelity, gradient and Hessian for a given waveform, specified in po- lar coordinates. Only the phase channel gradient is returned, the am- plitude profile is taken as a given. Syntax: [traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

## Parameters / inputs

- phi_profile -set of control pulse phases from an amplitude-phase
- description.

## Outputs

- fidelity -figure of merit for the overlap of the current state
- of the system and the desired state(s). When penalty
- methods are specified, fidelity is returned as an ar-
- ray separating the penalties from the simulation
- fidelity.
- gradient -gradient of the fidelity with respect to the control
- sequence. When penalty methods are specified, gradi-
- ent is returned as an array separating penalty gra-
- dients from the fidelity gradient.
- hessian -Hessian of the fidelity with respect to the control
- sequence. When penalty methods are specified, gradi-
- ent is returned as an array separating penalty Hes-
- sians from the fidelity Hessian.

## Implementation structure

- Cost function for optimal control using the GRAPE algorithm. Returns
- fidelity, gradient and Hessian for a given waveform, specified in po-
- lar coordinates. Only the phase channel gradient is returned, the am-
- plitude profile is taken as a given. Syntax:
- [traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)
- phi_profile -set of control pulse phases from an amplitude-phase
- description.
- fidelity -figure of merit for the overlap of the current state
- of the system and the desired state(s). When penalty
- methods are specified, fidelity is returned as an ar-
- ray separating the penalties from the simulation
- fidelity.
