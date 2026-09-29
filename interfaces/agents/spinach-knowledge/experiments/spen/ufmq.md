# experiments/spen/ufmq.m

Source: [experiments/spen/ufmq.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/ufmq.m)
Wiki: [ufmq.m](https://spindynamics.org/wiki/index.php?title=ufmq.m)
Reference: [DOI 10.1002/cphc.201800667](https://doi.org/10.1002/cphc.201800667), Figure 1A.

Signature: `fid=ufmq(spin_system,parameters,H,R,K,G,F)`

This routine implements the ultrafast multiple-quantum sequence in the `imaging()` context. That context supplies `H`, `R`, `K`, `G`, and `F`; the source requires the `sphten-liouv` formalism. The Liouvillian used for pulses and readout is `L=H+F+1i*R+1i*K`.

## State preparation and encoding

Starting from `parameters.rho0`, the sequence applies a 90-degree x pulse, evolves for `parameters.delay`, applies a 180-degree x pulse, and evolves for the same delay. The next 90-degree pulse is x for even `parameters.mqorder` and y for odd order. The code then selects coherence order +mqorder, applies a shaped chirp with positive encoding-gradient polarity, selects -mqorder, applies the second chirp with negative polarity, and selects +mqorder again. A third 90-degree x pulse converts to the selected +1 coherence for signal acquisition. These are the sequence operations in this implementation; no other pulse mechanism is implied.

The chirp waveform is generated from `pulsenpoints`, `Te`, `BW`, `nWURST`, and `chirptype`, with the source accepting `wurst` or `smoothed`. `Te` is the pulse duration in seconds, `BW` the bandwidth in Hz, and `Ge` the encoding gradient in T/m. The source divides `Te` into `pulsenpoints` equal time steps and uses the first encoding-gradient entry for the shaped pulses.

## Readout timing and output

The acquisition interval is `Taq=npoints*deltat`, with `deltat` in seconds and `Ga` the acquisition gradient in T/m. Before readout the routine prephases for `Taq/2` under `L-Ga*G{1}`. Each loop then evolves under the positive acquisition gradient for `Taq`, followed by `npoints` coil detections `coil'*rho` separated by `deltat` steps under the negative gradient. The state is carried from one loop to the next; the source uses nested `for` loops, not independent loop initial states. The returned complex `fid` has `npoints` rows and `nloops` columns.

The required parameter fields checked by the source are `rho0`, `coil`, `dims`, `npts`, `spins`, `npoints`, `nloops`, `Ga`, `deltat`, `pulsenpoints`, `BW`, `Ge`, `delay`, `Te`, `nWURST`, `chirptype`, and `mqorder`. `dims` is the sample dimension in metres; `npts` is the number of grid points; `npoints` is the number of acquired points per gradient readout; and `nloops` counts positive/negative readout loops. `spins` must be a one-entry cell array, which names the active nucleus. The source documents `Ga` and `Ge` in T/m, `deltat` and `Te` in seconds, and `BW` in Hz. The header also documents `offset` in Hz. When invoked through `imaging()`, the context applies this offset to `H` via `frqoffset` before `ufmq` receives the generator, so it is effective once even though `ufmq` does not read the field directly; do not add it again by hand.

The pulse train and FID shape above describe the source-defined simulation design; no experimental measurement or numerical run is reported here. The source header's synopsis spells the function name `ufmq_nmr`, while its MATLAB declaration is `ufmq`; the signature above follows the declaration.
