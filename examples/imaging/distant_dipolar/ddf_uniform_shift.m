% Uniform proton magnetisation after a 60-degree x pulse. With no
% gradients, diffusion, or relaxation, only the shape-dependent mean
% field survives. The fixed L+ receiver measures signed phase drift
% relative to the first sample in the resonant rotating frame.

function ddf_uniform_shift()

% Spin system in a 7 T field
sys.magnet=7;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};

% Complete spherical tensor basis
bas.formalism='sphten-liouv';
bas.approximation={'none'};
sys.parallel={'processes',1};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Homogeneous periodic cell with no flow, diffusion, or relaxation
parameters.dims=0.01;
parameters.npts=16;
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
parameters.npoints=21;
parameters.ddf_steps=4;
time_axis=(0:parameters.npoints-1)'/parameters.sweep;

% Compare signed phase drifts in three limiting sample shapes
parameters.flip_angles=[pi/3 0];
parameters.gradient_areas=[0 0];
factors=[1/3 0 1];
phase=zeros(parameters.npoints,3);
for k=1:3

    % Subtract the initial receiver phase before estimating the slope
    parameters.demag_factor=factors(k);
    fid=imaging(spin_system,@crazed_ddf,parameters);
    phase(:,k)=unwrap(angle(fid/fid(1)));

end

% Report the signed angular frequency and plot the accumulated phase
disp(table(factors.',(time_axis\phase).','VariableNames',{'Nz','rad_per_s'}));
kfigure();
plot(time_axis,phase);
kxlabel('time, s');
kylabel('phase relative to first sample, rad');
klegend('sphere','axial needle','transverse disc');

end


