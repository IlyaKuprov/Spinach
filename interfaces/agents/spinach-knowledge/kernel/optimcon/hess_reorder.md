# kernel/optimcon/hess_reorder.m

- Signature: `hess=hess_reorder(hess,K,N)`
- Source: [kernel/optimcon/hess_reorder.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/hess_reorder.m)

## Purpose

Reorders both variable axes of a Hessian to match the alternate flattening of a control waveform. If a waveform has K control rows and N time points, the input ordering is channel-fast within each time point, for example `[X1 Y1 Z1 X2 Y2 Z2 ... Xn Yn Zn]`. The output ordering is time-fast within each channel, for example `[X1 X2 ... Xn Y1 Y2 ... Yn Z1 Z2 ... Zn]`. Here the letter identifies a channel and the index identifies a time point. The same permutation converts in the reverse direction when supplied dimensions are correspondingly swapped; the routine does not infer the ordering.

## Inputs and output

- `hess` must be a numeric square matrix of size `(K*N) × (K*N)`.
- `K` and `N` must each be a positive integer scalar: respectively the number of waveform control rows and time points.
- The output `hess` is the permuted matrix. The gradient is not changed by this function; the source notes that its dimensions and element order follow the waveform.

The implementation reshapes `hess` to `[K N K N]`, permutes the axes as `[2 1 4 3]`, and reshapes to `[N*K N*K]`. The checks do not require `hess` to be real or symmetric.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=hess_reorder.m)
