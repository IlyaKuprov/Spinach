% A standard diffusion equation solver with no spin dynamics 
% present. Non-uniform diffusion coefficient distribution.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk

function diffusion_test_2b()

% Load the phantom
load('phantom_a.mat','R1');

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
parameters.dims=[0.02 0.02];
parameters.npts=[108 90];
parameters.deriv={'period',7};

% 2D flow parameters
parameters.u=zeros(parameters.npts);
parameters.v=zeros(parameters.npts);

% 2D diffusion tensor field
parameters.dxx=5e-5*kron(ones(parameters.npts(1),1),...
                         linspace(0,1,parameters.npts(2)).^2);
parameters.dxy=zeros(parameters.npts);
parameters.dyx=zeros(parameters.npts);
parameters.dyy=5e-5*kron(ones(parameters.npts(1),1),...
                         linspace(0,1,parameters.npts(2)).^2);

% Diffusion and flow generator
F=v2fplanck(spin_system,parameters); F=inflate(F);

% Timing parameters
timestep=5e-4; nsteps=200;

% System trajectory
traj=evolution(spin_system,F,[],R1(:),timestep,nsteps,'trajectory');

% Plotting
figure();
for n=1:nsteps
    imagesc(reshape(traj(:,n),108,90));
    drawnow; pause(0.025);
end

end

