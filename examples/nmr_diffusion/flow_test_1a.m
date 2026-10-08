% A standard diffusion and flow equation solver with no spin
% dynamics present and periodic boundary condition. Diffusion
% coefficient is constant throughout the sample.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk

function flow_test_1a()

% Ghost spin
sys.magnet=0;
sys.isotopes={'G'};

% No spin interactions
inter.zeeman.matrix=cell(1);
inter.coupling.matrix=cell(1);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sample geometry
parameters.dims=0.02;
parameters.npts=100;
parameters.deriv={'period',7};

% Diffusion and flow parameters
parameters.u=0.3*ones(parameters.npts,1);
parameters.dxx=5e-5*ones(parameters.npts,1);

% Diffusion and flow generator
F=v2fplanck(spin_system,parameters); F=inflate(F);

% Initial condition
rho=exp(-0.125*((1:100)-20).^2)';

% Timing parameters
timestep=5e-4; nsteps=70;

% System trajectory
traj=evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory');

% Physically correct X axis
x_axis=linspace(-parameters.dims/2,parameters.dims/2,parameters.npts);

% Plotting
kfigure();
for n=1:nsteps
    plot(x_axis',traj(:,n)); kgrid;
    kxlabel('sample coordinate, m'); kylabel('concentration');
    axis([-parameters.dims/2 parameters.dims/2 0 1]);
    drawnow; pause(0.025);
end

end

