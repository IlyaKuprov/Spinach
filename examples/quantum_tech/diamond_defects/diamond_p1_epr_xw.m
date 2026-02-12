% Field-swept EPR spectra of P1 center in diamond at X and W bands. Syntax:
%
%                        diamond_p1_epr_xw()
%
% This function takes no input arguments.
%
% This function produces diagnostic plots.
% alexey.bogdanov@weizmann.ac.il

function diamond_p1_epr_xw()

% Check consistency.
grumble();

% Set P1 model parameters.
p1_params.temperature=296;
p1_params.concentration=0.001;
p1_params.orientation='111';

% Build the spin system.
[sys,inter]=diamond_p1(p1_params);
%[sys,inter]=diamond_p1_15n(p1_params);

% Set the magnetic field for fieldsweep.
sys.magnet=1;

% Define the basis set.
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Run Spinach housekeeping.
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set common EPR parameters.
parameters.spins={'E'};
parameters.grid=6;
parameters.fwhm=1e-4;
parameters.int_tol=0.01;
parameters.tm_tol=0.1;
parameters.npoints=1024;
parameters.rspt_order=Inf;

% Set X-band parameters.
parameters.mw_freq=9.5e9;
parameters.window=[0.33 0.35];

% Run the X-band simulation.
[b_axis_x,spec_x]=fieldsweep(spin_system,parameters);

% Compute the X-band derivative spectrum.
spec_x=field_deriv(b_axis_x,spec_x);

% Plot the X-band spectrum.
kfigure(); plot(b_axis_x',spec_x');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('P1 X-band EPR');
axis tight; kgrid;

% Set W-band parameters.
parameters.mw_freq=94e9;
parameters.window=[3.348 3.36];

% Run the W-band simulation.
[b_axis_w,spec_w]=fieldsweep(spin_system,parameters);

% Compute the W-band derivative spectrum.
spec_w=field_deriv(b_axis_w,spec_w);

% Plot the W-band spectrum.
kfigure(); plot(b_axis_w',spec_w');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('P1 W-band EPR');
axis tight; kgrid;

end

% Consistency enforcement.
function grumble()

% This function has no inputs to validate.

end

% Finite-difference derivative on an evenly-spaced grid.
function spec_der=field_deriv(b_axis,spec)

% Set the axis spacing.
delta_b=b_axis(2)-b_axis(1);

% Preallocate the output.
spec_der=zeros(size(spec),'like',1i);

% Apply central differences.
spec_der(2:end-1)=(spec(3:end)-spec(1:end-2))/(2*delta_b);

% Apply one-sided differences at the ends.
spec_der(1)=(spec(2)-spec(1))/delta_b;
spec_der(end)=(spec(end)-spec(end-1))/delta_b;

end
