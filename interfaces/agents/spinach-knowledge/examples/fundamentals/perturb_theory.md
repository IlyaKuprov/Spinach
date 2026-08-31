# examples/fundamentals/perturb_theory.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/perturb_theory.m`
- Signature: `perturb_theory()`
- Total lines: 81

## Purpose

Rayleigh-Schrodinger and Van Vleck perturbation theory modules test. Eigenvector representations differ in the two theories, but the energies are the same.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Settings; implemented by `ham_dim=512`.
- Lines 14-15: H0 -Zeeman interaction; implemented by `sigma=pauli(ham_dim)`.
- Lines 18-19: H1 -random matrix; implemented by `H1=randn(ham_dim)+1i*randn(ham_dim)`.
- Lines 22-23: Energies -perturbation theories; implemented by `E_rs=zeros(ham_dim,max_ord)`.
- Lines 32-33: Energies -diagonalisation; implemented by `E_inf=eig(H0+H1,'vector')`.
- Lines 36-37: Comparison with diagonalisation; implemented by `diffs_rs=[diag(H0) real(E_rs)]-E_inf`.
- Lines 48-49: Eigensystems, PTs; implemented by `V_rs=cell(max_ord,1)`.
- Lines 53-54: RSPT gets the eigensystem directly; implemented by `[~,V_rs{n}]=rspert(diag(H0),H1,n)`.
- Lines 56-57: VVPT returns a generator that needs exponentiation; implemented by `[~,V_vv{n}]=vvpert(diag(H0),H1,n); V_vv{n}=expm(V_vv{n})`.
- Lines 61-62: Zero order is unit matrices; implemented by `V_rs=[{eye(ham_dim)}; V_rs]`.
- Lines 65-66: Eigensystem -diagonalisation; implemented by `[V_inf,E_inf]=eig(H0+H1,'vector')`.
- Lines 70-71: Comparison with diagonalisation; implemented by `diffs_rs=cellfun(@(V)norm(abs(V'*V_inf)-eye(size(V)),2),V_rs)`.

### Control flow inferred from the code

- Line 25: `for` loop over `n=1:max_ord`.
- Line 51: `for` loop over `n=1:max_ord`.

### Key state/data transformations

- Lines 10: computes `ham_dim` using `ham_dim=512`.
- Lines 11: computes `max_ord` using `max_ord=10`.
- Lines 12: computes `per_amp` using `per_amp=1/25`.
- Lines 15: computes `sigma` using `sigma=pauli(ham_dim)`.
- Lines 16: computes `H0` using `H0=full(sigma.z)`.
- Lines 19: computes `H1` using `H1=randn(ham_dim)+1i*randn(ham_dim)`.
- Lines 23: computes `E_rs` using `E_rs=zeros(ham_dim,max_ord)`.
- Lines 24: computes `E_vv` using `E_vv=zeros(ham_dim,max_ord)`.
- Lines 26: computes `E_rs(:,n)` using `E_rs(:,n)=rspert(diag(H0),H1,n)`.
- Lines 28: computes `E_vv(:,n)` using `E_vv(:,n)=vvpert(diag(H0),H1,n)`.
- Lines 33: computes `E_inf` using `E_inf=eig(H0+H1,'vector')`.
- Lines 37: computes `diffs_rs` using `diffs_rs=[diag(H0) real(E_rs)]-E_inf`.
- Lines 39: computes `diffs_vv` using `diffs_vv=[diag(H0) real(E_vv)]-E_inf`.
- Lines 49: computes `V_rs` using `V_rs=cell(max_ord,1)`.
- Lines 50: computes `V_vv` using `V_vv=cell(max_ord,1)`.
- Lines 54: computes `[~,V_rs{n}]` using `[~,V_rs{n}]=rspert(diag(H0),H1,n)`.
- Lines 57: computes `[~,V_vv{n}]` using `[~,V_vv{n}]=vvpert(diag(H0),H1,n); V_vv{n}=expm(V_vv{n})`.
- Lines 66: computes `[V_inf,E_inf]` using `[V_inf,E_inf]=eig(H0+H1,'vector')`.

## Implementation structure

- Rayleigh-Schrodinger and Van Vleck perturbation theory
- modules test. Eigenvector representations differ in the
- two theories, but the energies are the same.
- Settings
- H0 -Zeeman interaction
- H1 -random matrix
- Energies -perturbation theories
- Energies -diagonalisation
- Comparison with diagonalisation
- Eigensystems, PTs
- RSPT gets the eigensystem directly
- VVPT returns a generator that needs exponentiation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pauli()`, `E_rs()`, `rspert()`, `E_vv()`, `vvpert()`, `kfigure()`, `subplot()`, `set()`, `kxlabel()`, `kylabel()`, `klegend()`, `V_inf()`, `cellfun()`.
