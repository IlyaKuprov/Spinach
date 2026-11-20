% Uniform-field phase and an analytic gradient echo in a 1D sample.
% Equal voxel weights give the discrete counterpart of a sinc echo.
%
% talos@spindynamics.org

function field_profile_echo()

% Set up one proton and a uniform 20 Hz resonance offset
sys.isotopes={'1H'};
sys.magnet=5.9;
sys.parallel={'processes',1};
sys.output='hush';
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
spin_system=basis(create(sys,inter),bas);

% Set a stationary, equally weighted sample and one gradient reversal
parameters.spins={'1H'};
parameters.offset=20;
parameters.dims=0.04;
parameters.npts=41;
parameters.deriv={'period',3};
parameters.g_step_dur=1e-4;
parameters.g_n_steps=20;
parameters.rlx_ph={};
parameters.rlx_op={};
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(parameters.npts,1)/parameters.npts};
parameters.coil_st={coil_state(spin_system,'L+','1H')};
parameters.u=0;
parameters.diff=0;

% Establish the uniform-field phase, then apply the linear field profile
parameters.g_amp=0;
uniform=imaging(spin_system,@grad_echo,parameters);
parameters.g_amp=1e-4;
echo=imaging(spin_system,@grad_echo,parameters);

% Calculate discrete analytic phase and refocusing at the actual voxel centres
duration=parameters.g_n_steps*parameters.g_step_dur;
times=(-parameters.g_n_steps:parameters.g_n_steps)*parameters.g_step_dur;
positions=linspace(-parameters.dims/2,parameters.dims/2,parameters.npts);
profile=mean(exp(1i*spin('1H')*parameters.g_amp*positions.'*times),1);
phase=exp(-2i*pi*parameters.offset*(times+duration));
assert(max(abs(uniform(:)/uniform(1)-phase(:)))<1e-7);
assert(max(abs(echo(:)./uniform(:)-profile(:)))<1e-7);

% Plot the complex signal after removing the uniform-field phase
kfigure();
plot(times,real(echo(:)./uniform(:)),times,real(profile),'--');
kxlabel('time from echo centre, seconds');
kylabel('normalised in-phase signal');
klegend({'Spatial calculation','Analytic voxel sum'});

end


