# examples/fundamentals/spin_lock.m

- Signature: `spin_lock()`
- Source: [`examples/fundamentals/spin_lock.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/spin_lock.m)

## Model and question

This two-proton example follows the Cartesian magnetisation components during a spin-lock period. It illustrates coherent evolution of the coupled pair under a transverse locking field; it is not an orientation or powder-quadrature calculation, and the source specifies no state-space truncation test.

The system is two `1H` spins at 5.9 T, with scalar shifts 1.0 and 1.5 and a 7.0 Hz scalar coupling. It uses the full `sphten-liouv` basis (`bas.approximation='none'`). The initial operator is `4*state(...,'Lz','all')`; an x-axis `pi/2` rotation prepares the transverse state. The Hamiltonian is the NMR-assumption Hamiltonian plus the y-directed spin-lock term `2*pi*1.5e3*Ly`. No relaxation term is added.

## Propagation and output

`evolution` uses multichannel mode with a 1e-4 s time step and 100 steps (nominally 0.01 s). Its six observables are `Lx`, `Ly`, and `Lz` for spin 1 followed by the same components for spin 2. The example plots the real-valued component triplets as two trajectories on a Bloch-sphere surface.

The source defines no pass/fail tolerance, convergence criterion, or comparison target; the sphere limits are plot bounds, not acceptance thresholds. The code displays trajectories but does not report a numerical fit or claim a measured spin-lock efficiency.
