# examples/quantum_tech/circuit_qed/cavity_fock_grape_b.m

- Signature: `cavity_fock_grape_b()`
- Source: [cavity_fock_grape_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/cavity_fock_grape_b.m)

## Purpose

This example uses band-limited GRAPE to transfer a cavity from vacuum to its two-photon Fock level while leaving a dispersively coupled qubit in its upper level. It optimises smooth control waveforms rather than independent piecewise-constant control amplitudes.

## Model and control objective

The system is a four-level truncated cavity mode (`C4`) and a spin-1/2 electron (`E`) at zero magnet field. The cavity is on resonance with the control frame (`inter.modes.frqs={0 []}`), and the declared cavity–qubit dispersive coupling is 656.2 kHz. Spinach converts the mode-frequency and coupling inputs from Hz to angular-frequency units internally; pulse times are seconds, and the control scale is applied directly to the Hamiltonian in angular-frequency units.

The cavity-QED assumption uses a common rotating frame and rotating-wave approximation, retaining the full dispersive term χ N Lz: cavity photon number couples to the qubit’s Lz component and shifts its transition frequency. The four controls are the two cavity quadratures and the two qubit rotations (`Lx`, `Ly`). The normalised initial state uses cavity level `BL1` (vacuum) and qubit level `ZL2`; the target changes the cavity to `BL3` (two photons) while retaining `ZL2`. No dissipative generator is added, so this is a closed-system coherent-state-transfer calculation, not a cavity-preparation measurement.

## Numerical objective and reported output

Each of the 40 slices lasts 33 ns, giving a total pulse duration of 1.32 μs. The two sine and three cosine basis functions span the 40-point waveform; the four control channels are expanded in this common basis. The code sets the control scale to 1.76828×10^7 rad/s, uses the `NS` penalty with weight 0.001, and calls limited-memory BFGS through `fmaxnewton` for at most 300 iterations. The source header describes runtime on a minutes scale; this is not a timing result from this review.

After optimisation, the script reconstructs the waveform and propagates all slices directly. Its transfer score is `real(trace(rho_targ' * rho))`; the script errors if this score is below 0.80. The plot converts the control scale to MHz by dividing by 2π×10^6. The source file contains no recorded score or convergence history.

The source header attributes its model and parameters to the matching example in the paraqeet package.

## Scope

The model is the example’s truncated, closed cavity–qubit system. It demonstrates a simulated Fock-state transfer under smooth controls; it is not evidence of experimental Fock-state preparation, device performance, or convergence beyond the source’s stated acceptance check.
