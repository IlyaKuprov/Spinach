% Field-swept powder EPR spectra of an NV centre 
% in diamond at X and W bands.
%
% alexey.bogdanov@weizmann.ac.il

function diamond_nvm_epr_xw()

% Set NV centre parameters
nv_params.temperature=296;
nv_params.concentration=0.001;
nv_params.orientation='111';
nv_params.nitrogen='14N';

% Build the spin system
[sys,inter]=diamond_nvm_gs(nv_params);

% Field sweep
sys.magnet=1;

% Define the basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% EPR sim parameters
parameters.spins={'E3'};
parameters.grid=20;
parameters.fwhm=0.001;
parameters.int_tol=0.1;
parameters.tm_tol=0.01;
parameters.npoints=512;
parameters.rspt_order=Inf;

% Set X-band parameters
parameters.mw_freq=9.5e9;
parameters.window=[0.1 0.5];

% Run the X-band simulation
[b_axis_x,spec_x]=fieldsweep(spin_system,parameters);

% Plot the X-band spectrum
kfigure(); scale_figure([1.50 0.75]);
subplot(1,2,1); plot(b_axis_x',spec_x');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('NV$^{-}$ X-band EPR');
xlim tight; ylim padded; kgrid;

% Set W-band parameters
parameters.mw_freq=94e9;
parameters.window=[3.2 3.5];

% Run the W-band simulation
[b_axis_w,spec_w]=fieldsweep(spin_system,parameters);

% Plot the W-band spectrum
subplot(1,2,2); plot(b_axis_w',spec_w');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('NV$^{-}$ W-band EPR');
xlim tight; ylim padded; kgrid;

end

