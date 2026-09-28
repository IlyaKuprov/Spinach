# kernel/overloads/@ttclass/ranks.m

- Signature: `ttranks=ranks(ttrain)`

## Purpose

Returns the bond dimensions of the tensor trains stored in the input buffer.

## Physical / mathematical content

For each buffered train, the function reads each core's first dimension as the corresponding left rank and the final core's fourth dimension as the right boundary rank.

## Numerical / algorithmic content

The output is an `(ncores+1)-by-ntrains` array; the first and last entries for each train are 1 for valid tensor trains.

## Parameters / inputs

- ttrain - a tensor train object

## Outputs

- ttranks - `(ncores+1)` by `ntrains` array; the first and last elements for each train are 1

## Implementation structure

- Read the number of cores and buffered trains.
- Allocate the output array.
- For each train, extract left ranks from the cores and the final right rank from the last core.
