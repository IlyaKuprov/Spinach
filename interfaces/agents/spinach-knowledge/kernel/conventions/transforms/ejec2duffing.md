# kernel/conventions/transforms/ejec2duffing.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/ejec2duffing.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ejec2duffing.m) · [Koch et al.](https://doi.org/10.1103/PhysRevA.76.042319)

## Conversion

Using the asymptotic transmon expressions cited by the source, ejec2duffing maps Josephson energy ej and charging energy ec to the Duffing oscillator frequency and anharmonicity:

~~~text
frq    = sqrt(8*ej.*ec) - ec
anharm = -ec
~~~

Both inputs are energy divided by Planck's constant, in Hz. frq is in Hz for inter.modes.frqs; anharm is in Hz for inter.modes.anharms. The expressions are applied elementwise, so both outputs have the common input size. The source states that the approximation is accurate only deep in the transmon regime ej/ec >> 1; a warning is issued if any element has ej/ec < 20.

## Inputs and constraints

ej and ec must be numeric, real, finite, positive arrays with identical sizes. No particular dimensionality or nonempty-array requirement is checked. The ratio threshold is a warning condition, not a rejection condition.
