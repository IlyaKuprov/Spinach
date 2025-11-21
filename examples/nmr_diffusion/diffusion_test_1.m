% A standard diffusion equation solver with no spin
% dynamics present.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk

function diffusion_test_1()

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
parameters.u=zeros(parameters.npts,1);
parameters.diff=5e-5;

% Diffusion and flow generator
F=v2fplanck(spin_system,parameters); F=inflate(F);

% Initial condition
rho=exp(-0.125*((1:100)-20).^2)';

% Timing parameters
timestep=5e-4; nsteps=90;

% System trajectory
traj=evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory');

% Physically correct X axis
x_axis=linspace(-parameters.dims/2,parameters.dims/2,parameters.npts);

% Plotting
kfigure(); 
for n=1:nsteps
    plot(x_axis',traj(:,n)); kylabel('concentration'); 
    kxlabel('sample coordinate, m'); kgrid;
    axis([-parameters.dims/2 parameters.dims/2 0 1]);
    drawnow; pause(0.025);
end

end

