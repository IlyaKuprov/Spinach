# kernel/kinetics/kinetics.m

- Signature: `K=kinetics(spin_system)`

## Purpose

Builds the chemical-reaction and magnetization-flux kinetics superoperator from the chemical settings in `spin_system.chem`.

## Physical / mathematical content

Reaction rates transfer basis-state populations between source and destination spin systems. Optional magnetization fluxes transfer single-spin orders and handle multi-spin orders according to the selected intramolecular or intermolecular model. Optional radical-pair recombination adds singlet/triplet kinetics using the selected `rp_theory` model.

## Numerical / algorithmic content

The routine assembles sparse transitions from the configured reaction rates, checks compatibility of the source and destination basis subspaces, and adds the configured flux and radical-pair terms. Chemical reactions and magnetization-flux handling require the `sphten-liouv` formalism; radical-pair recombination accepts `sphten-liouv` or `zeeman-liouv`.

## Parameters / inputs

- `spin_system` - Spinach spin-system description, including the chemical reaction rates, magnetization flux settings, and (when enabled) radical-pair theory and rates.

## Outputs

- `K` - chemical kinetics superoperator. When assembling a Liouvillian manually, include it as `1i*K`, for example `L=H+1i*R+1i*K`. Spinach context functions include kinetics automatically.

## Implementation structure

The routine initializes a sparse superoperator, adds configured chemical-reaction transitions, processes magnetization fluxes when present, and then adds the selected radical-pair recombination model when enabled.
