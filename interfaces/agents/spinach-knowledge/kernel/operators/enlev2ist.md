# kernel/operators/enlev2ist.m

- Signature: `[states,coeffs]=enlev2ist(mult,lvl_num,particle)`

## Purpose

Expands the projector onto a specified Zeeman energy level as a linear combination of irreducible spherical tensors.

## Physical / mathematical content

The energy-level numbering depends on the particle type: spin levels are counted from the bottom up, while bosonic levels are counted from the top down. The resulting projector is represented in the irreducible spherical tensor basis.

## Numerical / algorithmic content

The routine builds a diagonal projector of size mult by mult. For a spin it sets element (mult-lvl_num+1,mult-lvl_num+1) to one; for a boson it sets (lvl_num,lvl_num) to one. It then calls oper2ist to obtain the tensor-basis states and coefficients.

## Parameters / inputs

- mult - multiplicity of the spin or dimension of the bosonic level space; a positive integer.
- lvl_num - energy-level number, counted from the bottom for spins and from the top for bosons, within the range 1 to mult.
- particle - particle type: 'S' for a spin or 'B' for a boson.

## Outputs

- states - indices in the Spinach IST basis that contribute to the operator; use lin2lm to convert them to L,M spherical-tensor indices.
- coeffs - coefficients of the irreducible spherical tensors in the linear combination.

## Implementation structure

1. Check the multiplicity and level-number inputs.
2. Place a unit entry on the appropriate diagonal position for the selected particle type; other particle values produce an error.
3. Expand the projector with oper2ist and return its states and coefficients.
