# experiments/overtone/overtone_a.m

MATLAB source: [experiments/overtone/overtone_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/overtone/overtone_a.m)

This function acquires an overtone signal in the frequency domain. Its input `spins` is a single-element cell array naming the overtone-active nucleus. The function does not implement overtone cross-polarisation or a pseudocontact-tensor calculation; it prepares an overtone frequency and delegates acquisition to `slowpass`.

## Inputs and frequency mapping

`spectrum=overtone_a(spin_system,parameters,H,R,K)` requires `spins`, `sweep`, `npoints`, `rho0`, and `coil`, together with the context-supplied matrices `H`, `R`, and `K`. `sweep` is a two-element frequency interval in Hz around the overtone frequency, and `npoints` is the requested positive integer sample count.

The source computes `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`, replaces the sweep with `ovt_frq-parameters.sweep`, then calls `slowpass`. Thus the frequency remapping is explicit in the code; the mapped endpoints are passed to the acquisition routine without additional pulse or tensor construction in this function.

## Output and limits

`spectrum` is the frequency-domain signal for the supplied starting state and detection coil, with `npoints` samples. The source notes that relaxation must be present for the `slowpass` matrix inversion and that `R` should not be thermalised. It does not specify a field-orientation or MAS sweep here.

https://spindynamics.org/wiki/index.php?title=overtone_a.m
