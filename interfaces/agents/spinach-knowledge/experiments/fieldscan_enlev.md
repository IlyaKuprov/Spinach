# experiments/fieldscan_enlev.m

- Signature: `fieldscan_enlev(spin_system,parameters)`

## Purpose

Plots a user-specified number of the lowest energy levels of the system as a function of the applied magnetic field. The energies are obtained using the Arnoldi method. Syntax: fieldscan_enlev(spin_system,parameters)

## Physical / mathematical content


- At a fixed Euler-angle orientation, constructs the field-dependent Hamiltonian from the oriented Zeeman term and coupling contribution.
- Computes and plots the real energy levels versus magnetic field; energy values are converted to cm^-1 for display.

## Numerical / algorithmic content


- Samples a linear magnetic-field grid spanning `parameters.fields` and computes the Hamiltonian at each point.
- The optional `eigs` calculation runs only when `parameters.nstates` is supplied, selecting the requested states for an active-space projection; no time propagation is performed.

## Parameters / inputs

- parameters.fields -two-element vector in Tesla,
- ordered as [from to]
- parameters.npoints -number of points in the scan
- parameters.orientation -system orientation, three-
- element vector containing
- Euler angles in radians,
- ordered as [alp bet gam]
- parameters.nstates -number of lowest energy states
- to solve for

## Implementation structure


- Validates the `zeeman-hilb` setup, field range, grid size, orientation, and optional state count, then constructs the oriented Zeeman and coupling terms.
- Evaluates the requested energies across the field grid and plots them against field.
