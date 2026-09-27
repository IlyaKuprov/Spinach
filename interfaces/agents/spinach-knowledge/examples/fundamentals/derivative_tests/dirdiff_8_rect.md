# examples/fundamentals/derivative_tests/dirdiff_8_rect.m

- Signature: `dirdiff_8_rect()`

## Purpose

Tests the phase-modulated GRAPE Hessian by finite-differencing analytical gradients for a rectangular pulse.

## Setup

The test covers sphten-liouv, zeeman-liouv, and zeeman-hilb, using a 13C system from dirdiff_test_system. The channel map is [1;1;1;1]. The four control operators are Lx, Ly, 0.4*Lx+0.2*Ly, and 0.7*Ly-0.1*Lx; initial states are Sx, Sy, Sz, with targets −Sz, Sy, Sx. Power levels are 2*pi*linspace(50e3,70e3,10). GRAPE uses Newton optimisation, a 1000-iteration limit, and the rectangle integrator. Five intervals have pulse_dt=12.8e-6*ones(1,5); the amplitude rows are ones(1,5) and 0.8+0.1*(1:5).

## Hessian check

For a random 2-by-5 phase array randn(2,5)/3, the code requests the analytical Hessian from grape_phase. It uses centered differences of gradients with h=1e-5, perturbing waveform entries i=1, i=end, and i=5 to test the leftmost, rightmost, and fifth Hessian columns. Each selected column passes when norm(hess_anl(:,i)-hess_num,1)<1e-5*norm(hess_num,1); otherwise the test raises an error naming the formalism and column.
