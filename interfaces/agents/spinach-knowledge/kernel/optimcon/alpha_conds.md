# kernel/optimcon/alpha_conds.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/alpha_conds.m
Wiki: https://spindynamics.org/wiki/index.php?title=alpha_conds.m

## Purpose and interface

`alpha_conds(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)` evaluates one selected line-search acceptance condition and returns its logical result. It tests a candidate; it does not perform the search or propagate a state.

## Conditions

- `test_type=0`: strict monotonic increase, `fx_1>fx_0`.
- `test_type=1`: Armijo sufficient increase, `fx_1>=fx_0+spin_system.control.ls_c1*alpha*(gfx_0'*dir)`.
- `test_type=2`: strong Wolfe curvature, `abs(gfx_1'*dir)<=spin_system.control.ls_c2*abs(gfx_0'*dir)`.
- `test_type=3`: positive trial directional derivative, `gfx_1'*dir>0`.

The comparisons and constants above are exactly those in the branch expressions. The routine reads `ls_c1` for type 1 and `ls_c2` for type 2; the gradient-direction products use MATLAB transpose, with the relevant vectors checked as real columns.

## Input checks

`test_type` must be a real numeric scalar in `0:3`. Types 0 and 1 require real scalar `fx_0` and `fx_1`; type 1 additionally requires real scalar `alpha`, real column `gfx_0`, and real column `dir`. Type 2 requires real columns `gfx_0`, `gfx_1`, and `dir`; type 3 requires real columns `gfx_1` and `dir`. The source checks matching dimensions between the direction and the gradients used in each condition. `spin_system` supplies the control coefficient only for types 1 and 2.
