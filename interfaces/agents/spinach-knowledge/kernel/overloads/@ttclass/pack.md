# kernel/overloads/@ttclass/pack.m

- Signature: `ttout=pack(tt)`
- Source: [`kernel/overloads/@ttclass/pack.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/pack.m)
- Wiki: [`ttclass/pack.m`](https://spindynamics.org/wiki/index.php?title=ttclass/pack.m)

## Core action

For `N==1` buffered train, returns the input unchanged. Otherwise, it allocates new zero-filled cores with physical dimensions from `tt.sizes` and internal bond dimensions equal to the sums of corresponding ranks across buffered trains; the two boundary ranks are fixed to one. For trains with more than one site, it places the coefficient-scaled first core for each summand in its own first-bond block, places intermediate cores on rank-diagonal blocks, and places each last core in its own final-bond block. This block/direct-sum construction represents the sum without expanding it to a dense matrix. For the one-site case, it instead adds each coefficient-scaled core into the sole output core.

The packed train is materialised in its core arrays and is not recompressed: the summed internal ranks remain, even if a smaller representation could exist. Coefficients are multiplied into the first core (or the one-site core) without conjugation.

## Output and guards

For multiple buffered trains, the output has `coeff=1`, the combined cores, and `tolerance=abs(sum(tt.tolerance))`; the logical site sizes are unchanged. A single buffered train is returned as-is, including its existing coefficient and tolerance. The implementation has no explicit type, shape, or rank validation. The source recommends using `ttclass/shrink.m` rather than calling `pack` directly in normal use.
