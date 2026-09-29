# kernel/utilities/rspert.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rspert.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rspert.m)

## Purpose

`rspert.m` implements Rayleigh–Schrödinger perturbation theory to arbitrary order for a non-degenerate Hamiltonian `H0 + H1`, returning perturbative corrections to eigenvalues and eigenvectors. The implementation follows Eqs. 2.21–2.23 from Stefan Stoll's PhD thesis.

## Behaviour

- Syntax: `[Ep,Vp]=rspert(E0,H1,order)`.
- A consistency-checking subfunction `grumble` validates the inputs (real column vector `E0`, Hermitian `H1`, consistent dimensions, positive integer `order`).
- Reciprocal energy differences `Q = 1./(E0'-E0)` are computed, with the diagonal zeroed; if any element of `Q` is non-finite, the function errors with `H0 has degenerate energy levels.`
- First order: `E{1}=diag(H1)` and `V{1}=Q.*H1`.
- Higher orders (loop `k=2:order`): computes `R=H1*V{k-1}`, sets `E{k}=real(diag(R))`, and builds `V{k}` by subtracting lower-order products `V{k-m}.*E{m}'` for `m=1:(k-1)` before multiplying elementwise by `Q`.
- Summation: `Ep` starts from `E0` and accumulates all `E{n}` for `n=1:order`; `Vp` starts from the identity and accumulates all `V{n}`.
- Normalisation: `Vp` is column-normalised as `Vp./sqrt(sum(abs(Vp).^2,1))`.
- Notes from the header: there must be no degeneracies in `H0`; `H1` must be Hermitian; the source header cautions that perturbation theory requires `norm(H1,2)` much smaller than the smallest energy gap in `H0`, that numerical artefacts can appear beyond sixth order, and that the stated complexity is linear in order and cubic in matrix dimension.

## Inputs and outputs

**Inputs**

- `E0` — eigenvalues of `H0`, a column vector of real numbers.
- `H1` — perturbation, written in the basis that diagonalises `H0`.
- `order` — order of perturbation theory to be used; 6 is the sensible maximum.

**Outputs**

- `Ep` — eigenvalues of `H0+H1` to the specified order, a vector of reals, not necessarily sorted in the same way as the input.
- `Vp` — normalised eigenvectors of `H0+H1` to the specified order in perturbation theory, a square unitary matrix with eigenvectors in columns, in the same order as the eigenvalues in `Ep`.

## References

- Stefan Stoll's PhD thesis, Eqs. 2.21–2.23.
- Spinach Wiki page: [rspert.m](https://spindynamics.org/wiki/index.php?title=rspert.m)
