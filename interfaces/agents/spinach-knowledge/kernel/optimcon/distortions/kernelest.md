# kernel/optimcon/distortions/kernelest.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/kernelest.m`
- Signature: `h=kernelest(x,y,ker_len,method,align,lambda)`
- Total lines: 134

## Purpose

FIR convolution kernel estimation from input and output signal samples. Syntax: h=kernelest(x,y,ker_len,method,align,lambda)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- x -input samples on a uniform grid
- y -output samples on the same grid
- ker_len -kernel length (number of taps)
- method -'backslash' (default) | 'pinv' | 'svd' | 'tikh'
- align -'causal' (default) or 'same' output alignment
- lambda -Tikhonov parameter for 'tikh' (optional)

## Outputs

- h -estimated convolution kernel

## Implementation structure

- FIR convolution kernel estimation from input and output signal
- samples. Syntax:
- h=kernelest(x,y,ker_len,method,align,lambda)
- x -input samples on a uniform grid
- y -output samples on the same grid
- ker_len -kernel length (number of taps)
- method -'backslash' (default) | 'pinv' | 'svd' | 'tikh'
- align -'causal' (default) or 'same' output alignment
- lambda -Tikhonov parameter for 'tikh' (optional)
- h -estimated convolution kernel
- Set the defaults
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `toeplitz()`, `lower()`, `conv_mat()`, `eps()`, `s_inv()`, `s_vals()`, `isvector()`, `isscalar()`, `ischar()`, `isstring()`.
