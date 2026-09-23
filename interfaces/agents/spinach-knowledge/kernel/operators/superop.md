# kernel/operators/superop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/superop.m`
- Signature: `A=superop(spin_system,opspec,side)`
- Total lines: 215

## Purpose

Sided product superoperator in the spherical tensor basis set. Returns superoperators corresponding to right or left multiplication of a density matrix by a user-specified operator. Syntax: A=superop(spin_system,opspec,side)

## Physical / mathematical content

- The product of a basis operator with the density matrix expanded in irreducible spherical tensor products is read off the Lie structure tables of the single-spin algebras (`spin_system.bas.lpst` and `rpst`, one page per operator index of each multiplicity): the operator specified by `opspec` maps each source state of the active spins to destination states with structure coefficients, and the multi-spin table is the direct product of the single-spin tables.
- Commutation superoperators (`leftofcomm`, `rightofcomm`) drop the paths in which either the source or the destination is the unit state on every active spin, because those terms cancel between the left and right products; `comm` and `acomm` are assembled from the two sided calls.

## Numerical / algorithmic content

- For every source pattern of the active spins, the basis states carrying that pattern are located by column comparisons on the sparse descriptor `spin_system.bas.basis`, and the same is done for the destination pattern; when the two subspaces coincide state for state, the superoperator elements are written directly, otherwise the source rows are matched to the destination rows with `ismember(...,'rows')` over the columns of the passive spins.
- The row matching runs on a dense copy of the passive columns held in the smallest signed integer class that covers every state index of the system (`min_int_type` of the largest multiplicity squared minus one), built once per call from the non-zeros of the sparse descriptor; integer rows sort and match several times faster than sparse double rows, and the column comparisons that select the subspaces stay on the sparse descriptor.
- The result is returned in XYZ triplet form (row, column, value) for the caller to assemble; an empty operator is returned as a single zero triplet.
- The grumbler requires a `sphten-liouv` basis and an integer opspec row with one entry per spin, each below the squared multiplicity of that spin.

## Parameters / inputs

- opspec - Spinach operator specification described in Sections 2.1 and 3.3 of http://dx.doi.org/10.1016/j.jmr.2010.11.008
- side - 'left' or 'right' causes the function to return a product superoperator corresponding to a product from that side; 'comm' or 'acomm' results in commutation and anticommutation superoperator respectively

## Outputs

- A - a three-column array of row indices (first column), column indices (second column) and values (third column)
- Note: this is a very general function to which direct calls are not usually required, please use the (much friendlier) operator() function.
- Note: the superoperator is returned in XYZ sparse format, which is different from Matlab's CSC format.
