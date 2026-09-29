# kernel/conventions/transforms/ppm2hz.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/ppm2hz.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ppm2hz.m)

- Signature: `hz=ppm2hz(ppm,B0,nucleus)`

## Behaviour

The conversion is `hz = 1e-6 * ppm * (B0 * spin(nucleus) / (2*pi))`. The factor `1e-6` converts ppm to a dimensionless fraction; `spin(nucleus)` supplies the signed magnetogyric ratio in rad/(s*T), and division by `2*pi` converts angular frequency to cycles per second (Hz). The sign is preserved.

## Inputs

- `ppm` — real numeric scalar or array of chemical shifts in ppm.
- `B0` — real numeric scalar magnetic induction in tesla.
- `nucleus` — character array naming an isotope; the source gives `'1H'` as an example.

## Output

- `hz` — resonance offset in Hz, with the same array shape as `ppm` (or a scalar when `ppm` is scalar).

The explicit checks require numeric real `ppm`, numeric real scalar `B0`, and a character-array `nucleus`. The function imposes no additional shape restriction on `ppm`.
