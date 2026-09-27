# examples/optimal_control/features_curv.m

- Signature: `features_curv()`

## Purpose

Demonstrates optimal-control pulse design in curvilinear coordinates. The control variables are mapped to RF amplitude and phase through a user-defined coordinate map and its Jacobian, which are passed to the curvilinear GRAPE calculation.

## Physical / mathematical content

The example uses a two-13C spin system at 14.1 T with a 60 Hz scalar coupling and the sphten-liouv basis. It explores the control landscape through curvilinear coordinates for the RF amplitude and phase.

## Numerical / algorithmic content

The script supplies the coordinate transformation and Jacobian to grape_curv and optimises over an ensemble of 11 B1 power levels spanning 0.6 to 1.4. This is a gradient-based curvilinear-control calculation. For the coherence-to-singlet transfer, the source gives the optimal fidelity as 1/√2 ≈ 0.7071.

## Implementation structure

The MATLAB code constructs the spin system and basis, defines the amplitude/phase coordinate map and its Jacobian, sets the B1 power ensemble, and calls grape_curv for pulse optimisation. The ensemble uses linspace(0.6,1.4,11).
