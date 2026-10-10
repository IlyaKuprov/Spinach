# kernel/optimcon/distortions/amp_root.m

Source: [kernel/optimcon/distortions/amp_root.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/amp_root.m)
Wiki: [Spinach reference](https://spindynamics.org/wiki/index.php?title=amp_root.m)

## Purpose

Applies a smooth radial amplitude compression independently to every X/Y control pair in every time slice. For a pair <code>(X,Y)</code>, let <code>r = sqrt(X^2+Y^2)</code>, saturation level <code>a</code>, and sharpness <code>s</code>. Both components are multiplied by <code>[1+(r/a)^s]^(-1/s)</code>, so the output radial amplitude is <code>r/[1+(r/a)^s]^(1/s)</code>. The pair's direction is unchanged. A starting choice of <code>s=4</code> is documented.

The waveform uses rad/s nutation-frequency units. The saturation level has the same amplitude units and sets the limiting radial output amplitude. The routine is a distortion map, not an objective or constraint function.

## Call and data

<code>[w,J] = amp_root(w,sat_lvls,s)</code>

Rows of <code>w</code> are successive X and Y components, and each column is one time slice. <code>sat_lvls</code> and <code>s</code> each provide one value per X/Y pair. There are no default input values. <code>w</code> must be real numeric with an even row count; each saturation level must be finite, real and positive; each sharpness value must be a real positive integer. The implementation checks element counts against the number of pairs but does not explicitly require finite entries in <code>w</code> or finite <code>s</code> values.

The second output <code>J</code> is optional. When requested, it is a sparse <code>numel(w)</code>-by-<code>numel(w)</code> Jacobian with respect to MATLAB column-major vectorisation of <code>w</code>. Each pair contributes a 2-by-2 block of the form <code>scale*I + curvature*[X;Y]*[X Y]</code>, where <code>scale = [1+(r/a)^s]^(-1/s)</code> and <code>curvature = -r^(s-2)/a^s * [1+(r/a)^s]^(-1/s-1)</code>. At <code>r=0</code>, the implementation explicitly sets <code>scale=1</code> and <code>curvature=0</code>. If needed, Jacobian block values are gathered to host memory before sparse assembly.
