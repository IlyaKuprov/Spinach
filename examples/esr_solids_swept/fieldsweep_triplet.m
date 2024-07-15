% Powder averaged X-band field-swept ESR spectrum of photo-
% generated pentacene triplet state.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk
% guinevere.mathies@uni-konstanz.de

function fieldsweep_triplet()

% Magnet field (must be 1)
sys.magnet=1;

% Triplet electron
sys.isotopes={'E3'};

% Zeeman tensor, assumed isotropic
inter.zeeman.matrix={diag([2.0 2.0 2.0])};

% ZFS, photo-excited pentacene triplet
D=1360.1*1e6; E=-47.2*1e6;    % [Hz]
inter.coupling.matrix={zfs2mat(D,E,0,0,0)};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E3'};
parameters.grid=4;
parameters.mw_freq=9e9;           % Hz
parameters.fwhm=2e-3;             % Tesla
parameters.int_tol=0.01;
parameters.tm_tol=0.01;
parameters.window=[0.25 0.40];  % Tesla
parameters.npoints=512;
parameters.rspt_order=Inf;

% Zeeman tensor into Hz/Tesla
Z=-spin('E')*inter.zeeman.matrix{1,1}/(2*pi*2.00231930436256);

% Orientation- and field-dependent intial condition
parameters.rho0=@(B,alp,bet,gam)zftrip(spin_system,euler2dcm(alp,bet,gam)*...
                                                   inter.coupling.matrix{1,1}*... 
                                                   euler2dcm(alp,bet,gam)',...
                                                   [0.56 0.31 0.13],...
                                                   euler2dcm(alp,bet,gam)*Z*...
                                                   euler2dcm(alp,bet,gam)',B,1);

% Run the simulation
[b_axis,spec]=fieldsweep(spin_system,parameters);

% Plotting
figure(); plot(b_axis',spec');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.'); kgrid;
xlim tight; drawnow;

end

