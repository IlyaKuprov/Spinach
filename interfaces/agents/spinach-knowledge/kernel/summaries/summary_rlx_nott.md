# kernel/summaries/summary_rlx_nott.m

- Signature: `summary_rlx_nott(spin_system)`

## Behaviour

Reports four stored Nottingham DNP relaxation rates: electron R1 and R2, then nuclear R1 and R2, from `spin_system.rlx.nott_r1e`, `nott_r2e`, `nott_r1n` and `nott_r2n`. The report labels each value in Hz; this routine prints the values directly without converting them.

The function returns nothing and routes its lines through `report(spin_system,...)`. Its local guard requires `spin_system` to be a structure.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_rlx_nott.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_rlx_nott.m)
