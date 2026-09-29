# examples/singlet_states/singlet_imaging_1.m

Source: [examples/singlet_states/singlet_imaging_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/singlet_imaging_1.m)

## Purpose

Compare singlet-order and transverse-magnetisation images in a simulated sample with one-dimensional flow and diffusion. The example first prepares transverse magnetisation with a shaped RF pulse, converts it to singlet order with an M2S block, then propagates both states in the imaging model.

## Spin system and imaging model

The system is a pair of 13C spins at 9.4 T, with scalar offsets 0.03 and -0.03 and a 55 Hz scalar coupling. Their coordinates are entered as [0, 0, 0] and [1.2, 0, 0]. Redfield relaxation is configured with zero equilibrium, secular retention, and a 1 ns correlation time; both relaxation tolerances are 1e-5. The sample grid has dimensions [0.10, 0.015] and 150 by 15 points. The spatial derivative option is period 7, and the sequence offset is 0. A tube-shaped phantom occupies columns 6 through 10; the relaxation phantom is uniform.

The imaging Liouvillian combines the Hamiltonian, flow, relaxation, and diffusion terms. The flow field is u = -0.06 and v = 0; the diffusion tensor has only an x component, 3.6e-6. The source does not attach units to these spatial-field entries. A 6 mT/m x-gradient is included during RF excitation. The displayed spatial extent is set in the plotting code and labelled in millimetres.

## RF preparation and singlet conversion

The excitation uses a Gaussian amplitude shape with 25 steps over 2 ms, a 2 kHz pulse frequency, and phase pi/2. Each step lasts 80 microseconds; the amplitude table is 2*pi*100 times the normalised Gaussian shape. shaped_pulse_af applies the pulse with the x-gradient term. The resulting transverse-magnetisation state is saved before conversion.

The M2S block uses J = 55 and delta-v = 6 to set its free-evolution interval and number of repeated blocks: m1 is floor(pi*J/(2*delta_v)), incremented by one if odd. Alternating free evolution under the Hamiltonian with x rotations produces the singlet-order state used for the second image. The source also propagates the saved transverse-magnetisation state as a comparison.

## Propagation and observable

Each state is propagated for 100 steps of 0.01 s in the imaging context. The singlet channel is detected with the two-spin singlet operator; the magnetisation channel uses the total 13C L+ state. The script reconstructs each phantom at every trajectory step and displays the real singlet-order and transverse-magnetisation images side by side. The plotted singlet range is 25*[-2e-3, 6e-3]/2 and the transverse-magnetisation range is 25*[-2e-3, 6e-3], labelled in millimetres. This describes the requested simulation and plotting operations, not a validated image or measured lifetime.

## Citation

The source file does not identify a DOI or publication citation.
