# examples/spin_chemistry/singlet_yield_4.m

- Signature: `singlet_yield_4()`

## Purpose

Reproduce Figure 3 from the paper by Till, Timmel, Brocklehurst and Hore: http://dx.doi.org/10.1016/S0009-2614(98)01158-0. The source notes that the original paper uses only electron Zeeman operators for the field sweep, missing effects associated with the increasing nuclear Zeeman interaction on the high-field side of the plot. Calculation time: seconds.

## Physical / mathematical content

- The spin system contains two electrons (`E`) and two protons (`1H`). Electron scalar Zeeman values are `2.0023` and `2.0044`; both proton scalar Zeeman entries are zero.
- Scalar couplings between electron 1 and protons 3 and 4 are `gauss2mhz(35)*1e6` and `gauss2mhz(30)*1e6`, respectively. The `{4,4}` scalar coupling entry is set to zero.
- The calculation plots singlet recombination yield against `log(magnetic induction / mT)`.

## Numerical / algorithmic content

- Set `sys.magnet=1` for the field sweep. Use the `zeeman-hilb` formalism with approximation `none`.
- Sweep `parameters.fields=1e-3*10.^linspace(-5,3,2000)` with rates `[0.1 1.0 10.0 100.0 1000.0]*1e6`. Set `parameters.electrons=[1 2]`, `parameters.spins={'E'}`, and `parameters.needs={'zeeman_op'}`.
- Create the spin system, apply the basis, and run `liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Define the unit magnet, spin system, basis, couplings, and sequence parameters.
- Initialise the Spinach spin system with `create` and `basis`.
- Simulate with `liquid`, then plot the result against `linspace(-5,3,2000)` with grid and axis labels.
