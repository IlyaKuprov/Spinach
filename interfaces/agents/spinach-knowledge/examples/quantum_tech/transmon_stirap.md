# examples/quantum_tech/transmon_stirap.m

Source: [examples/quantum_tech/transmon_stirap.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/transmon_stirap.m)

## What it models

An ensemble-robust GRAPE pulse optimisation for transfer from BL1 to BL3 of a single three-level Duffing transmon. The Gaussian STIRAP-like pair is the initial pulse guess; the source estimates calculation time in minutes, runs an optimisation, and does not report a measured transfer or establish agreement with an experiment.

## Model and pulse guess

The source cites [10.1038/ncomms10628](https://doi.org/10.1038/ncomms10628) for the model and parameters. It sets the field to zero, uses a T3 mode with rotating-frame ladder detuning 10e6 Hz and anharmonicity -20e6 Hz, and builds the cavity/Duffing drift. The basis is unapproximated Zeeman-Liouville. The two 3-by-3 control matrices couple levels 0-1 and 1-2 through their respective ladder quadratures. Their commutator superoperators encode coherent Hamiltonian control in Liouville space; no relaxation superoperator is added.

The initial guess samples t from -150 to 150 ns at 300 points. Both pulses have Gaussian width parameter sigma=45 ns; the 1-2 pulse is centred at ts/2=-45 ns and precedes the 0-1 pulse centred at 0, with ts=-90 ns. The source sets the normalised pulse amplitudes from 43.4 MHz and 38.2 MHz and uses a power ensemble of 2*pi*[30,35,40,45,50]*1e6 rad/s. Each control sample has the source-defined adjacent-time spacing, about 1.003 ns.

## Optimisation and scope

The initial and target states are BL1 and BL3, respectively, each normalised by its Frobenius norm. GRAPE is configured for the two controls with the SNSA penalty, weight 1.0, Goodwin method, and a 50-iteration limit; fmaxnewton is called with the Gaussian pair as its initial pulse. Optimisation plotting is enabled for controls, spectrogram, and robustness. These settings define an optimisation task, not a reported final fidelity. No relaxation, rotational diffusion, correlation spectrum, cross-correlations, or secular settings are present; the Liouville formalism alone does not imply relaxation.
