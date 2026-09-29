# kernel/summaries/summary_rlx_weiz.m

- Signature: `summary_rlx_weiz(spin_system)`

## Behaviour

Prints the stored Weizmann DNP relaxation rates `weiz_r1e`, `weiz_r2e`, `weiz_r1n` and `weiz_r2n` with electron/nuclear R1/R2 labels, plus each nonzero entry of the inter-nuclear dipolar rate matrices `weiz_r1d` and `weiz_r2d` with its row and column indices. The report labels the rates in Hz and performs no conversion.

There is no return value: output goes through `report(spin_system,...)`. The local guard requires `spin_system` to be a structure.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_rlx_weiz.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_rlx_weiz.m)
