% Second-order CRAZED echo of a synthetic scalar-coupled proton pair.
% Two distinct shifts evolve with isotropic J coupling in the rotating
% frame. Feedback is from the total local magnetisation; the receiver
% averages the complex signal over one pitch. Both traces use the same
% resonant x-pulse calibration, with diffusion and relaxation absent.

function ddf_coupled_pair()

% Spin system in a 7 T field
sys.magnet=7;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0.02};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=7;

% Complete spherical tensor basis
bas.formalism='sphten-liouv';
bas.approximation={'none'};
sys.parallel={'processes',1};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Homogeneous periodic cell with no flow, diffusion, or relaxation
parameters.dims=0.01;
parameters.npts=32;
parameters.deriv={'period',3};
parameters.u=zeros(parameters.npts,1);
parameters.diff=0;
parameters.rlx_ph={};
parameters.rlx_op={};

% Longitudinal preparation and uniform volume-average receiver
parameters.spins={'1H'};
parameters.offset=0;
parameters.rho_ref=state(spin_system,'Lz','1H');
parameters.rho0=kron(ones(parameters.npts,1),parameters.rho_ref);
parameters.coil=kron(ones(parameters.npts,1)/parameters.npts,coil_state(spin_system,'L+','1H'));

% Set the equilibrium magnetisation by the chosen dipolar time in seconds
tau_d=0.1;
parameters.magnetisation=1/(spin('1H')*spin_system.tols.mu0*tau_d);
parameters.direction=[0 0 1];
parameters.sweep=200;
parameters.npoints=1;
parameters.ddf_steps=4;

% Calibrate the receiver with an isolated 90-degree x pulse
parameters.demag_factor=1;
parameters.flip_angles=[pi/2 0];
parameters.gradient_areas=[0 0];
ref=imaging(spin_system,@crazed_ddf,parameters);

% Wind one pitch and then two pitches across the cell
parameters.npoints=21;
time_axis=(0:parameters.npoints-1)'/parameters.sweep;
parameters.flip_angles=[pi/2 pi/2];
parameters.gradient_areas=[1 2]*2*pi/(spin('1H')*parameters.dims);
echo=imaging(spin_system,@crazed_ddf,parameters)/ref;

% Repeat with feedback disabled but the same molecular Hamiltonian
parameters.magnetisation=0;
control=imaging(spin_system,@crazed_ddf,parameters)/ref;

% Display both receiver quadratures on a common single-pulse scale
kfigure();
plot(time_axis,[real(echo) imag(echo) real(control) imag(control)]);
kxlabel('time after mixing, s');
kylabel('signal / single-pulse signal');
klegend('real, DDF on','imaginary, DDF on','real, DDF off','imaginary, DDF off');

end


