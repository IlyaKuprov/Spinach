# kernel/grids/grid_kron.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/grid_kron.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=grid_kron.m)

## Purpose

Form a direct-product spherical orientation grid by combining two input grids through active ZYZ Euler rotations.

## Inputs and outputs

- `angles1` and `angles2` are real matrices with three columns `[alpha beta gamma]`, in radians. Let their row counts be n1 and n2.
- `weights1` and `weights2` are finite real column vectors with n1 and n2 entries, respectively.
- The outputs are `angles`, with n1*n2 rows and three Euler-angle columns in radians, and `weights`, an n1*n2-entry column vector.

## Product rule

Each input orientation is converted to a quaternion. For each pair of rows (i,j), the output quaternion is the Hamilton product q1_i*q2_j, converted back to ZYZ Euler angles by `qter2euler`. The row order holds each row of grid 1 while iterating through all rows of grid 2. The corresponding weight is `weights1(i)*weights2(j)`; equivalently, the output vector is `kron(weights1,weights2)`.

The function does not normalise the input or output weights and imposes no positivity or unit-sum requirement. It checks that each angle input is a real numeric three-column matrix and that each weight input is a finite real column with one entry per angle row; it does not check angle finiteness or constrain Euler angles to a canonical range. Pairing and output order are deterministic.
