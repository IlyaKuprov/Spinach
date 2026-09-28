# kernel/conventions/transforms/weblab2nqi.m

- Signature: `varargout=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)`

Converts Weblab one-cone model parameters (see `weblab_cone.png`) into NQI quadrupolar coupling tensors used by Spinach. See the [function reference](https://spindynamics.org/wiki/index.php?title=weblab2nqi.m).

## Inputs

- `C_q`: quadrupolar coupling constant `e^2*q*Q/h`, in Hz.
- `eta_q`: quadrupolar tensor asymmetry parameter.
- `I`: spin quantum number; an integer or half-integer at least 1.
- `alpha`, `theta`, `phi`: Weblab cone-model angles in radians. `phi` is required for two- and three-output calls and must be omitted for four- and six-output calls.

## Output modes

Each output `Q1`, `Q2`, … is a 3×3 quadrupolar coupling tensor in Hz. The function calls `eeqq2nqi(C_q,eta_q,I,[azimuth theta alpha])` for each site:

| Outputs | Call arguments | Azimuths, in output order |
| --- | --- | --- |
| `Q1,Q2` | `C_q,eta_q,I,alpha,theta,phi` | `-phi/2`, `+phi/2` |
| `Q1,Q2,Q3` | `C_q,eta_q,I,alpha,theta,phi` | `-phi`, `0`, `+phi` |
| `Q1`–`Q4` | `C_q,eta_q,I,alpha,theta` | `0`, `pi/2`, `pi`, `3*pi/2` |
| `Q1`–`Q6` | `C_q,eta_q,I,alpha,theta` | `0`, `pi/3`, `2*pi/3`, `pi`, `4*pi/3`, `5*pi/3` |

## Validation

Only 2, 3, 4, or 6 outputs are supported. All supplied inputs must be numeric, real scalars; `I` must also satisfy the spin constraint above. The function rejects a missing `phi` in two- or three-output mode and a supplied `phi` in four- or six-output mode.