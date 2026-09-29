# kernel/utilities/rspt_eig.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rspt_eig.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rspt_eig.m)

## Purpose

Computes the eigensystem of sparse Hamiltonians to a user-specified order in Rayleigh–Schrödinger perturbation theory (RSPT), with careful handling of diagonal dominance and an option to perform exact diagonalisation (expensive). The function also returns eigenvalue derivatives and transition moments between eigenvectors under a user-specified operator. The parametrisation matches use cases in field-swept EPR spectroscopy.

## Behaviour

- Syntax: `[E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)`.
- Input consistency is enforced by an internal `grumble` subfunction, which requires `parameters.rspt_order` to be present and to be a positive real scalar or `Inf`, requires `Hz`, `Hc`, and `Hmw` to be square numeric matrices, and requires `B` to be a real scalar.
- If the basis contains an `irrep` subfield (symmetry-adapted basis), the function loops over irreps: each Hamiltonian component (`Hz`, `Hc`) and the microwave operator `Hmw` are projected with the irrep projector `P` and symmetrised as `(X+X')/2`; a recursive call is made with the `irrep` field removed from the basis, and the resulting eigenvectors are projected back with `P`. Irrep blocks are concatenated.
- Otherwise, the Hamiltonian is formed as `H = B*Hz + Hc` and symmetrised (`full((H+H')/2)`), and the method is selected by `parameters.rspt_order`:
  - Orders 1–4: the Hamiltonian is split into diagonal and off-diagonal parts, and `rspert` is called with the specified order.
  - `Inf`: full diagonalisation via `eig(H,'vector')`.
  - Any other value: an error `'unsupported perturbation theory order.'` is raised.
- Energies are sorted in ascending order and eigenvectors reordered to match.
- If more than two output arguments are requested, `dE` is computed via the Hellmann–Feynman theorem as `real(diag(V'*Hz*V))`.
- If more than three output arguments are requested, transition moments are computed as `T = abs(V'*Hmw*V).^2`.
- If more than four output arguments are requested, level populations `LP = real(diag(V'*rho0*V))` are computed from the thermal equilibrium state:
  - If `parameters.rho0` is a function handle, it is called as `rho0(B, alpha, beta, gamma)` using `parameters.orientation(1:3)` (orientation-dependent equilibrium).
  - If `parameters.rho0` is a matrix, it is used directly (orientation-independent).
  - If `rho0` is not provided, the equilibrium is computed at the current temperature, orientation, and magnetic field via `equilibrium(spin_system,H)` with `H = B*Hz + Hc` symmetrised.

## Inputs and outputs

**Inputs:**

- `spin_system` — Spinach spin system object.
- `parameters` — Parameter struct:
  - `parameters.rspt_order` — perturbation theory order used to account for the off-diagonal part of the Hamiltonian; `Inf` for exact diagonalisation.
  - `parameters.rho0` — optional; a matrix specifying a user-defined thermal equilibrium state, or a function handle `f(B,alp,bet,gam)` returning the thermal equilibrium at each orientation and magnetic field; if absent, the thermal equilibrium is computed at the current temperature, orientation, and magnetic field.
- `Hz` — laboratory frame Hamiltonian containing only Zeeman terms at 1 Tesla.
- `Hc` — laboratory frame Hamiltonian containing all spin-spin couplings but no Zeeman terms.
- `Hmw` — observable operator without the amplitude prefactor (2-norm should be around 1).
- `B` — magnetic field, Tesla.

**Outputs:**

- `E` — column vector of energies, sorted in ascending order (rad/s).
- `V` — matrix with eigenvectors in columns, sorted left to right in the same order as the energies.
- `dE` — column vector of dE/dB derivatives, sorted in the same order as the energies.
- `T` — matrix of transition moments under `Hmw`.
- `LP` — column vector of energy level populations, sorted in the same order as the energies.

## References

- Spinach Wiki: [rspt_eig.m](https://spindynamics.org/wiki/index.php?title=rspt_eig.m)
- Source code: [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rspt_eig.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rspt_eig.m)
