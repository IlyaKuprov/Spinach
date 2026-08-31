# examples/dnp_liq/jdnp/energy_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/energy_levels.m`
- Signature: `energy_levels()`
- Total lines: 52

## Purpose

Energy level diagram transition from the Zeeman limit to the exchange coupling limit in a two-electron system.

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: 600 MHz magnet; implemented by `sys.magnet=14.1`.
- Lines 11-12: Two electrons; implemented by `sys.isotopes={'E','E'}`.
- Lines 14-15: Exaggerate g-factor difference; implemented by `inter.zeeman.scalar={1.9 2.1}`.
- Lines 17-18: Hilbert space calculation; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 21-22: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 25-26: Relevant operators; implemented by `Hz=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 31-32: omega_j scan from 0 to omega_e; implemented by `omega_e=sys.magnet*spin('E')`.
- Lines 35-36: Get the energy levels; implemented by `for n=1:numel(omega_j)`.
- Lines 38-39: Diagonalise the Hamiltonian; implemented by `[~,D]=eig(Hz-omega_j(n)*Lj,'vector')`.
- Lines 41-42: Sort and record energies; implemented by `energies(:,n)=sort(D)`.
- Lines 46-47: Plot energies; implemented by `kfigure(); plot(omega_j/omega_e,energies'/omega_e,'k-')`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=1:numel(omega_j)`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.9 2.1}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 26: computes `Hz` using `Hz=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 27-29: computes `Lj` using `Lj=operator(spin_system,{'Lx','Lx'},{1,2})+ operator(spin_system,{'Ly','Ly'},{1,2})+ operator(spin_system,{'Lz','Lz'},{1,2})`.
- Lines 32: computes `omega_e` using `omega_e=sys.magnet*spin('E')`.
- Lines 33: computes `omega_j` using `omega_j=linspace(-3*omega_e,3*omega_e,1000)`.
- Lines 39: computes `[~,D]` using `[~,D]=eig(Hz-omega_j(n)*Lj,'vector')`.
- Lines 42: computes `energies(:,n)` using `energies(:,n)=sort(D)`.

## Implementation structure

- Energy level diagram transition from the Zeeman limit to
- the exchange coupling limit in a two-electron system.
- 600 MHz magnet
- Two electrons
- Exaggerate g-factor difference
- Hilbert space calculation
- Spinach housekeeping
- Relevant operators
- omega_j scan from 0 to omega_e
- Get the energy levels
- Diagonalise the Hamiltonian
- Sort and record energies

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `spin()`, `omega_j()`, `energies()`, `kfigure()`, `kxlabel()`, `kylabel()`.
