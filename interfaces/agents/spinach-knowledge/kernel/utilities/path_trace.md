# kernel/utilities/path_trace.m

- Signature: `projectors=path_trace(spin_system,L,rho)`

## Purpose

Liouvillian path tracing. Treats the user-supplied Liouvillian as the adjacency matrix of a graph, computes the weakly connect- ed subgraphs of that graph and returns a cell array of project- ors into independently evolving populated subspaces. Syntax: projectors=path_trace(spin_system,L,rho)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- L -Hamiltonian or Liouvillian matrix
- rho -the initial state (source state screening)
- or the detection state (destination state
- screening); pass [] to disable screening

## Outputs

- projectors -a cell array of projectors into independently
- evolving populated subspaces. The projectors
- are to be used as follows:
- L_reduced=P'*L*P; (for matrices)
- rho_reduced=P'*rho; (for state vectors)
- Note: further information on how this function works is availa-
- ble in our papers on this subject

## Implementation structure

- Liouvillian path tracing. Treats the user-supplied Liouvillian
- as the adjacency matrix of a graph, computes the weakly connect-
- ed subgraphs of that graph and returns a cell array of project-
- ors into independently evolving populated subspaces. Syntax:
- projectors=path_trace(spin_system,L,rho)
- L - Hamiltonian or Liouvillian matrix
- rho - the initial state (source state screening)
- or the detection state (destination state
- screening); pass [] to disable screening
- projectors -a cell array of projectors into independently
- evolving populated subspaces. The projectors
- are to be used as follows:
