# cce_capabilities — registered experiment contracts

`contract=cce_capabilities(pulse_sequence)` is the native registry consulted by `clexpn` before retained generators are built. It compares actual function handles; a matching name or five-argument signature does not register a callback.

| Experiment | Projection | Dissipation | Reduction and outputs |
|---|---|---|---|
| `free_evolution` | Generalized or projected | Yes | Multiplicative signal/observable; generalized observed density |
| `hahn_echo` | Generalized or projected | Yes | Multiplicative signal/observable; generalized observed density |
| `finite_drive` | Generalized only | Yes | Multiplicative signal/observable; observed density |
| `cpmg_echo` | Generalized or projected | No | Multiplicative signal/observable |
| `bath_pumped_hahn` | Generalized or projected | No | Multiplicative signal/observable |
| `cce_eseem` | Generalized or projected | No | Multiplicative signal/observable |
| `noise_correlation` | Generalized only | Bath-only jumps | Additive correlation |

Registered representations are `zeeman-hilb` and `zeeman-liouv`, subject to each experiment and generator restriction. `cce_generators` requires retained Liouville space for native generalized jumps and separately constructs conditional projected dissipative generators.

The registry returns the experiment handle, evaluator, supported formalisms, projection/dissipation flags, reduction algebra, density-reduction flag, outputs and contract version. The evaluator constructs the retained subsystem and generators and calls the five-argument experiment. It produces internal raw records for `cce_schedule` and `cce_reduce`; these records are not the public first output of `clexpn`.

For multiplicative experiments the raw record carries signal, observed-density and detector data where applicable, absolute error bounds, validity and propagation diagnostics. Hahn records include the actual physical retained IDs and raw response evidence. Noise records retain value, moment, mean_time, mean_zero, valid, and abs_error_bound for additive ensemble reduction.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_capabilities.m)

The noise evaluator is bound to noise_correlation, even when called directly without noise parameters. A foreign schedule or an explicitly different output is rejected with clexpn:noise before subsystem construction. Its contract version is distinct from the general dynamics evaluator.

Free evolution and its finite-drive delegate use capability version 3: checkpoints from the earlier master-propagation arithmetic are not accepted as current raw records. The public experiment signatures are unchanged.

Direct registered evaluation checks retained time-history memory before constructing a subsystem or generator, including density components and retained density pages. Hahn capability version 4 covers native stock master-equation pulse propagation and generalized selective refocusing.
