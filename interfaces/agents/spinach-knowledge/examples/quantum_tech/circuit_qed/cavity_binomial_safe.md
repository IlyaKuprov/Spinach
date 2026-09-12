# examples/quantum_tech/circuit_qed/cavity_binomial_safe.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/cavity_binomial_safe.m`
- Signature: `cavity_binomial_safe()`
- Total lines: 232

## Purpose

Binomial bosonic code |0L>=(|0>+|4>)/sqrt(2), |1L>=|2> in a cavity dispersively coupled to a flux-tunable transmon ancilla, and the protection of its coherences from 1/f flux noise by a Stark-assis- ted flux-noise evasion (SAFE) drive on the transmon, Sec. 4.4.1 and Fig. 4.4(a,b) of Yunwei Lu's PhD thesis (Northwestern University, 2026). The flux noise dephasing rates of the code and error space coherences are comput

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `dressed_ens()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Transmon anharmonicity placing the sweet spot of Eq. (C.22) at -30 MHz for a 10 MHz drive, Hz; implemented by `anharm=-67e6`.
- Lines 24-25: Transmon-cavity exchange coupling and detuning giving a 0.5 MHz dispersive shift, Hz; implemented by `g_bc=86e6; delta_bc=1.414e9`.
- Lines 27-28: Dispersive shift and the sensitivities of the cavity terms to the transmon frequency; implemented by `chi=2*(g_bc/delta_bc)^2*anharm; sens_c=(g_bc/delta_bc)^2; sens_x=-4*anharm*g_bc^2/delta_bc^3`.
- Lines 30-31: Transmon frequency sensitivity to flux, rad/s per flux quantum; implemented by `dwb_dphi=2*pi*6e9`.
- Lines 33-34: Flux noise amplitude in flux quanta and the ultraviolet cutoff, Hz; implemented by `noise_amp=1e-5; f_uv=5e7`.
- Lines 36-37: Drive amplitude, rad/s, and the transmon-drive detunings to scan, Hz; implemented by `omega0=2*pi*10e6; detunings=-(20:1:80)*1e6`.
- Lines 39-40: Transmon and cavity relaxation times, seconds; implemented by `t1_b=50e-6; t1_c=20e-3`.
- Lines 42-43: Time step, number of steps, infidelity sampling stride, and trajectory count; implemented by `dt=1e-8; nsteps=30000; stride=300; ntraj=100`.
- Lines 45-46: Length of the noise synthesis grid and the infrared cutoff at its frequency resolution, Hz; implemented by `nlong=2^19; f_ir=1/(nlong*dt)`.
- Lines 48-49: Fock state pairs whose coherences matter for the code, Sec. 4.4.1; implemented by `pairs=[2 0; 4 2; 4 0; 3 0; 4 3; 2 1; 3 1]`.
- Lines 51-52: Three-level transmon and a five-level cavity, both in their own frames; implemented by `sys.magnet=0; sys.isotopes={'T3','C5'}`.
- Lines 60-61: Hilbert space twin for the dressed states; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 65-66: Hamiltonian without the detuning term, detuning, drive, and flux noise operators; implemented by `H_hilb=hamiltonian(assume(ss_hilb,'labframe'))`.
- Lines 71-72: Positions of the bare |g,n> states; implemented by `idx=zeros(1,5)`.
- Lines 77-78: Flux noise dephasing rates of the code coherences with and without the drive, Eq. (4.93) with |ln(omega_ir*t)|=4; implemented by `rates=zeros(size(pairs,1),numel(detunings)); dw=2*pi*1e3`.
- Lines 84-85: Operating detuning at the median of the rate minima, and the undriven rates there; implemented by `[~,min_idx]=min(rates,[],2); op_idx=round(median(min_idx)); delta_bd=detunings(op_idx)`.
- Lines 89-90: Report the suppression at the operating detuning; implemented by `disp(['rate minima at ' num2str(-detunings(min_idx)/1e6,'%.0f ') 'MHz, operating detuning ' num2str(-delta_bd/1e6) ' MHz'])`.
- Lines 93-94: Validate the sweet spot window; implemented by `if any(min_idx==1)||any(min_idx==numel(detunings))`.

### Control flow inferred from the code

- Line 73: `for` loop over `n=1:5`.
- Line 79: `for` loop over `n=1:numel(detunings)`.
- Line 94: conditional branch on `any(min_idx==1)||any(min_idx==numel(detunings))`.
- Line 97: conditional branch on `any(rates_off./rates(:,op_idx)<5)`.
- Line 120: `for` loop over `m=1:ntraj`.
- Line 138: `for` loop over `k=1:2`.
- Line 142: conditional branch on `drives(k)==0`.
- Line 151: `for` loop over `j=1:3`.
- Line 161: `for` loop over `n=1:ngrid`.
- Line 167: `for` loop over `m=1:ntraj`.
- Line 169: `for` loop over `n=1:nsteps`.
- Line 171: conditional branch on `mod(n,stride)==0`.
- Line 179: `for` loop over `n=1:numel(time_axis)`.
- Line 186: `for` loop over `t=1:3`.

### Key state/data transformations

- Lines 22: computes `anharm` using `anharm=-67e6`.
- Lines 25: computes `g_bc` using `g_bc=86e6; delta_bc=1.414e9`.
- Lines 28: computes `chi` using `chi=2*(g_bc/delta_bc)^2*anharm; sens_c=(g_bc/delta_bc)^2; sens_x=-4*anharm*g_bc^2/delta_bc^3`.
- Lines 31: computes `dwb_dphi` using `dwb_dphi=2*pi*6e9`.
- Lines 34: computes `noise_amp` using `noise_amp=1e-5; f_uv=5e7`.
- Lines 37: computes `omega0` using `omega0=2*pi*10e6; detunings=-(20:1:80)*1e6`.
- Lines 40: computes `t1_b` using `t1_b=50e-6; t1_c=20e-3`.
- Lines 43: computes `dt` using `dt=1e-8; nsteps=30000; stride=300; ntraj=100`.
- Lines 46: computes `nlong` using `nlong=2^19; f_ir=1/(nlong*dt)`.
- Lines 49: computes `pairs` using `pairs=[2 0; 4 2; 4 0; 3 0; 4 3; 2 1; 3 1]`.
- Lines 52: computes `sys.magnet` using `sys.magnet=0; sys.isotopes={'T3','C5'}`.
- Lines 53: computes `inter.modes.frqs` using `inter.modes.frqs={0 0}`.
- Lines 54: computes `inter.modes.anharms` using `inter.modes.anharms={anharm []}`.
- Lines 55: computes `inter.modes.lifetimes` using `inter.modes.lifetimes={t1_b t1_c}`.
- Lines 56: computes `inter.modes.kerr` using `inter.modes.kerr=cell(2,2); inter.modes.kerr{1,2}=chi`.
- Lines 57: computes `inter.temperature` using `inter.temperature=0`.
- Lines 58: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 61: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.

### Local helper functions

- Line 222: `dressed_ens()` — `function suscept=dressed_ens(H_hilb,noise_op,drive_op,omega0,idx,dw)`.
  - Representative operation: `energies=zeros(2,numel(idx)); shifts=[-dw dw]`.
  - Representative operation: `for m=1:2`.

## Implementation structure

- Binomial bosonic code |0L>=(|0>+|4>)/sqrt(2), |1L>=|2> in a cavity
- dispersively coupled to a flux-tunable transmon ancilla, and the
- protection of its coherences from 1/f flux noise by a Stark-assis-
- ted flux-noise evasion (SAFE) drive on the transmon, Sec. 4.4.1 and
- Fig. 4.4(a,b) of Yunwei Lu's PhD thesis (Northwestern University,
- 2026). The flux noise dephasing rates of the code and error space
- coherences are computed from the flux sensitivities of the dressed
- cavity transition frequencies as functions of the transmon-drive
- detuning, Eq. (4.93); at the common minimum the logical state |+L>
- is then propagated for 300 microseconds along 1/f flux noise tra-
- jectories under the Lindblad master equation, with and without the
- drive, and the decoherence-only infidelity of Eq. (4.94) and the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `idx()`, `state()`, `int2str()`, `dressed_ens()`, `detunings()`, `rates()`, `dressed()`, `pairs()`, `median()`, `num2str()`, `any()`.
