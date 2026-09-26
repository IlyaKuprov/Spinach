# examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m

- Signature: `kernel_from_transm()`

## Purpose

Response function extraction from the transmission profile of the HiPER instrument. We are sending a linear chirp via the AWG from 93.3 GHz to 94.7 GHz over 1400 ns and record- ing the power response (so a square root must be taken to obtain the amplitude). The two spikes show when the chirp starts and finishes. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Response function extraction from the transmission profile
- of the HiPER instrument. We are sending a linear chirp via
- the AWG from 93.3 GHz to 94.7 GHz over 1400 ns and record-
- ing the power response (so a square root must be taken to
- obtain the amplitude). The two spikes show when the chirp
- starts and finishes.
- Rob Hunter, Hassane el-Mkami, Graham Smith,
- Yujie Zhao, Shebha Anandhi Jegadeesan,
- Guinevere Mathies, Ilya Kuprov
- Load power data and convert to amplitude
- Plot as received
- Get apodisation weights
