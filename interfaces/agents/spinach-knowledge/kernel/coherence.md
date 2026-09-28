# kernel/coherence.m

- Signature: `rho=coherence(spin_system,rho,spec)`

## Purpose

Keeps only the specified coherence orders in a state vector or density matrix. This is useful as an analytical replacement for complicated phase cycles.

## Physical / mathematical content

- Coherence orders are determined from projection quantum numbers of the basis states. For each entry in `spec`, the function retains states whose summed coherence order on the selected spins matches one of the specified orders. All entries in `spec` must be satisfied.

## Numerical / algorithmic content

- The function builds a mask for each coherence specification, intersects the masks, zeros the excluded states, and restores the original shape of `rho`.

## Parameters / inputs

- `spin_system` — spin system specifying the basis, formalism, and spins.
- `rho` — a state vector or a horizontal stack thereof; in `zeeman-hilb`, a density matrix or a horizontal stack thereof.
- `spec` — a cell array specifying which coherence orders to keep on which spins. For example, `{{'13C',[1 -1]},{'1H',-1}}` keeps states with coherence order `((1 OR -1 on 13C) AND (-1 on 1H))`. Spins may be specified by isotope, spin number, `'electrons'`, `'nuclei'`, or `'all'`.

## Outputs

- `rho` — the state vector or density matrix with undesired coherence orders zeroed out.
- Note: this function requires `sphten-liouv`, `zeeman-liouv`, or `zeeman-hilb` formalism; Fokker-Planck direct products are supported in the Liouville space formalisms. In `zeeman-hilb`, density matrices are stretched into Liouville space, filtered there, and folded back.

## Implementation structure

- Checks the inputs and formalism, computes coherence orders for the basis, applies the intersected mask, and reshapes the result to the original dimensions. Reports a warning if the resulting magnetization is nearly zero.
