# examples/quantum_tech/circuit_qed/cavity_dual_rail_safe.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/cavity_dual_rail_safe.m`
- Signature: `cavity_dual_rail_safe()`
- Total lines: 101

## Purpose

Flux-noise dephasing rates of two dual-rail qubits whose flux-tuna- ble transmon-coupled rails share one transmon ancilla, and their suppression by a Stark-assisted flux-noise evasion (SAFE) drive on the transmon, Sec. 4.4.2 and Fig. 4.4(c) of Yunwei Lu's PhD thesis (Northwestern University, 2026). The logical dephasing rate of each dual-rail qubit is set by the flux sensitivity of the dressed sin- gle-photon transit

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Transmon anharmonicity, Hz; implemented by `anharm=-200e6`.
- Lines 23-24: Exchange couplings and detunings of the two transmon-coupled rails, Hz; implemented by `g_bc=[50e6 60e6]; delta_bc=[1.414e9 2.0e9]`.
- Lines 26-27: Dispersive shifts and the sensitivities of the rail terms to the transmon frequency; implemented by `chi=2*(g_bc./delta_bc).^2*anharm; sens_c=(g_bc./delta_bc).^2; sens_x=-4*anharm*g_bc.^2./delta_bc.^3`.
- Lines 29-30: Transmon frequency sensitivity to flux, rad/s per flux quantum, and the noise amplitude, flux quanta; implemented by `dwb_dphi=2*pi*6e9; noise_amp=1e-5`.
- Lines 32-33: Drive amplitude, rad/s, and the transmon-drive detunings to scan, Hz; implemented by `omega0=2*pi*10e6; detunings=-(20:0.5:50)*1e6`.
- Lines 35-36: Three-level transmon and the two transmon-coupled rails, all in their own frames; implemented by `sys.magnet=0; sys.isotopes={'T3','C3','C3'}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Hamiltonian without the detuning term, detuning, drive, and flux noise operators; implemented by `H_hilb=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 54-57: Positions of the bare |g,0,0>, |g,1,0>, and |g,0,1> states; implemented by `idx=[find(diag(state(spin_system,{'BL1','BL1','BL1'},{1,2,3}))>0.5) find(diag(state(spin_system,{'BL1','BL2','BL1'},{1,2,3}))>0.5) find(diag(state(spin_system,{'BL1','BL…`.
- Lines 59-60: Flux noise dephasing rates of the two rails under the drive, Eq. (4.93) with |ln(omega_ir*t)|=4; implemented by `rates=zeros(2,numel(detunings)); dw=2*pi*1e3; shifts=[-dw dw]`.
- Lines 73-74: Undriven rates from the bare dispersive sensitivities, Eq. (4.95); implemented by `rates_off=noise_amp*dwb_dphi*sqrt(2*4)*sens_c`.
- Lines 76-77: Locate the minima and the common operating point; implemented by `[~,min_idx]=min(rates,[],2); common=round(mean(min_idx))`.
- Lines 83-84: Validate the sweet spots; implemented by `if any(min_idx==1)||any(min_idx==numel(detunings))`.
- Lines 94-95: Plot the dephasing rates, Fig. 4.4(c); implemented by `kfigure(); semilogy(-detunings/1e6,1e-9*rates',-detunings([1 end])/1e6,1e-9*[1; 1]*rates_off,'--'); kgrid`.

### Control flow inferred from the code

- Line 61: `for` loop over `n=1:numel(detunings)`.
- Line 63: `for` loop over `m=1:2`.
- Line 65: `for` loop over `j=1:3`.
- Line 84: conditional branch on `any(min_idx==1)||any(min_idx==numel(detunings))`.
- Line 87: conditional branch on `abs(detunings(min_idx(1))-detunings(min_idx(2)))>5e6`.
- Line 90: conditional branch on `any(rates_off'./rates(:,common)<5)`.

### Key state/data transformations

- Lines 21: computes `anharm` using `anharm=-200e6`.
- Lines 24: computes `g_bc` using `g_bc=[50e6 60e6]; delta_bc=[1.414e9 2.0e9]`.
- Lines 27: computes `chi` using `chi=2*(g_bc./delta_bc).^2*anharm; sens_c=(g_bc./delta_bc).^2; sens_x=-4*anharm*g_bc.^2./delta_bc.^3`.
- Lines 30: computes `dwb_dphi` using `dwb_dphi=2*pi*6e9; noise_amp=1e-5`.
- Lines 33: computes `omega0` using `omega0=2*pi*10e6; detunings=-(20:0.5:50)*1e6`.
- Lines 36: computes `sys.magnet` using `sys.magnet=0; sys.isotopes={'T3','C3','C3'}`.
- Lines 37: computes `inter.modes.frqs` using `inter.modes.frqs={0 0 0}`.
- Lines 38: computes `inter.modes.anharms` using `inter.modes.anharms={anharm [] []}`.
- Lines 39: computes `inter.modes.kerr` using `inter.modes.kerr=cell(3,3); inter.modes.kerr{1,2}=chi(1); inter.modes.kerr{1,3}=chi(2)`.
- Lines 40: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 41: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `H_hilb` using `H_hilb=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 49: computes `num_b` using `num_b=operator(spin_system,'N',1)`.
- Lines 50: computes `drive_op` using `drive_op=operator(spin_system,'C',1)+operator(spin_system,'A',1)`.
- Lines 51-52: computes `noise_op` using `noise_op=num_b+sens_c(1)*operator(spin_system,'N',2)+sens_c(2)*operator(spin_system,'N',3)+ sens_x(1)*operator(spin_system,{'N','N'},{1,2})+sens_x(2)*operator(spin_syste…`.
- Lines 55-57: computes `idx` using `idx=[find(diag(state(spin_system,{'BL1','BL1','BL1'},{1,2,3}))>0.5) find(diag(state(spin_system,{'BL1','BL2','BL1'},{1,2,3}))>0.5) find(diag(state(spin_system,{'BL1','BL…`.
- Lines 60: computes `rates` using `rates=zeros(2,numel(detunings)); dw=2*pi*1e3; shifts=[-dw dw]`.

## Implementation structure

- Flux-noise dephasing rates of two dual-rail qubits whose flux-tuna-
- ble transmon-coupled rails share one transmon ancilla, and their
- suppression by a Stark-assisted flux-noise evasion (SAFE) drive on
- the transmon, Sec. 4.4.2 and Fig. 4.4(c) of Yunwei Lu's PhD thesis
- (Northwestern University, 2026). The logical dephasing rate of each
- dual-rail qubit is set by the flux sensitivity of the dressed sin-
- gle-photon transition frequency of its transmon-coupled rail, Eqs.
- (4.93) and (4.95); the rates are computed from the eigenvalues of
- the rotating frame Hamiltonian of Eq. (4.16) at a fixed drive amp-
- litude as functions of the transmon-drive detuning. Both rates go
- through a minimum in the same detuning window, so that one drive
- protects both qubits.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `chi()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `sens_c()`, `sens_x()`, `state()`, `detunings()`, `shifts()`, `vecs()`, `idx()`, `energies()`, `vals()`, `rates()`.
