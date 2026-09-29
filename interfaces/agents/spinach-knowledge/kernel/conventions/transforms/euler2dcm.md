# kernel/conventions/transforms/euler2dcm.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/euler2dcm.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=euler2dcm.m)

## Conversion and convention

The function converts Euler angles in radians to a 3x3 direction cosine matrix using the active ZYZ convention of Brink & Satchler, Fig. 1a. Each factor is a counterclockwise rotation about its indicated axis:

~~~text
Rz(theta) = [cos(theta)  -sin(theta)  0;  sin(theta)  cos(theta)  0;  0  0  1]
Ry(theta) = [cos(theta)  0  sin(theta);  0  1  0;  -sin(theta)  0  cos(theta)]
R = Rz(alpha)*Ry(beta)*Rz(gamma)
~~~

The documented applications are v_out = R*v_in for a 3x1 vector and A_out = R*A_in*R' for a 3x3 interaction tensor.

## Inputs and constraints

Call either R=euler2dcm([alpha beta gamma]) or R=euler2dcm(alpha,beta,gamma). In the one-input form the code takes the first three elements by linear indexing; it does not check that the input is a three-element vector, so later elements are unused and fewer than three elements fail during indexing. In the three-input form each angle must be a numeric real scalar. The source checks no angle range or finiteness. The output is a 3x3 direction cosine matrix.
