# kernel/overloads/@ttclass/pack.m

- Signature: `ttout=pack(tt)`

## Purpose

Packs the trains in the addition buffer into one tensor train without recompressing it. The source advises using `ttclass/shrink.m` rather than calling this function directly in normal use.

## Physical / mathematical content

The output combines the buffered summands into a single train whose bond dimensions are the sums of the corresponding ranks. The coefficients are incorporated into the first core; no recompression is performed.

## Numerical / algorithmic content

If there is only one buffered train, the function returns the input unchanged. Otherwise it allocates cores for the summed ranks and copies the buffered cores into the corresponding rank blocks.

## Parameters / inputs

- tt - tensor train object with unprocessed additions

## Outputs

- ttout - tensor train with the additions buffer absorbed into its cores, but not recompressed

## Implementation structure

- Read the train sizes and ranks; return immediately when there is one train.
- Sum the ranks across buffered trains and allocate the combined cores.
- Collect the buffered trains into the single output train.
