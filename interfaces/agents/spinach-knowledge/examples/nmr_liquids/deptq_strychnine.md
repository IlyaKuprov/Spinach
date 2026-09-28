# examples/nmr_liquids/deptq_strychnine.m

- Signature: `deptq_strychnine()`

## Purpose

DEPTQ135 experiment on strychnine. Calculation time: minutes

## Physical / mathematical content

- Simulates a one-dimensional DEPTQ135 experiment on the 1H/13C strychnine spin system. The signal is generated with Spinach’s liquid-state DEPTQ sequence and scalar-coupling-mediated transfer.
- Natural-abundance 13C isotopomers are simulated separately; the free-induction decay is exponentially apodised and Fourier transformed to give the plotted carbon spectrum.

## Numerical / algorithmic content

- Uses a sphten-liouv / IK-2 basis with scalar-coupling connectivity and proximity level 1; temperature is 298 K and the field is 5.9 T. Sequence parameters are sweep `10000`, offset `[5000 0]`, `npoints=2048`, `zerofill=8196`, `J=150`, and `beta=3*pi/4`.
- Iteration over 13C isotopomers is parallelised with `parfor`; no GPU execution is present.

## Implementation structure

- Create the 1H/13C strychnine spin system; set the 5.9 T field and 298 K temperature, then configure the scalar-coupling basis and DEPTQ135 parameters.
- Generate 13C isotopomers and simulate each in parallel. Exponentially apodise (`exp`, 6), Fourier transform, and plot the real 13C spectrum.
