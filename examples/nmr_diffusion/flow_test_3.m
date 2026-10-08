% Circular flow in three-dimensional space in the absence
% of spin dynamics.
%
% Calculation time: minutes, faster on GPU.
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk

function flow_test_3()

% Ghost spin
sys.isotopes={'G'};
sys.magnet=5.9;

% No spin interactions
inter.zeeman.matrix=cell(1);
inter.coupling.matrix=cell(1);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.enable={'greedy'}; % 'gpu'

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sample geometry
parameters.dims=[0.02 0.02 0.02];
parameters.npts=[50 50 50];
parameters.deriv={'period',7};

% Get a 3D grid
[X,Y,Z]=ndgrid(linspace(-parameters.dims(1)/2,parameters.dims(1)/2,parameters.npts(1)),...
               linspace(-parameters.dims(2)/2,parameters.dims(2)/2,parameters.npts(2)),...
               linspace(-parameters.dims(3)/2,parameters.dims(3)/2,parameters.npts(3)));

% Get circular wind vectors
parameters.u=-1000*Y;
parameters.v=+1000*X;
parameters.w=zeros(size(X));

% Constant diffusion tensor field
parameters.dxx=8e-6*ones(parameters.npts);
parameters.dxy=zeros(parameters.npts);
parameters.dxz=zeros(parameters.npts);
parameters.dyx=zeros(parameters.npts);
parameters.dyy=8e-6*ones(parameters.npts);
parameters.dyz=zeros(parameters.npts);
parameters.dzx=zeros(parameters.npts);
parameters.dzy=zeros(parameters.npts);
parameters.dzz=8e-6*ones(parameters.npts);

% Get the initial condition
X0=[-0.003  0.001  0.003]; 
Y0=[-0.003  0.001 -0.003]; 
Z0=[-0.003  0.001  0.003]; 
sigma=2e-6; A=zeros(size(X));
A=A+(1/sqrt((2*pi)^3*sigma^3))*exp(-((X-X0(1)).^2+(Y-Y0(1)).^2+(Z-Z0(1)).^2)/(2*sigma));
A=A-(1/sqrt((2*pi)^3*sigma^3))*exp(-((X-X0(2)).^2+(Y-Y0(2)).^2+(Z-Z0(2)).^2)/(2*sigma));
A=A+(1/sqrt((2*pi)^3*sigma^3))*exp(-((X-X0(3)).^2+(Y-Y0(3)).^2+(Z-Z0(3)).^2)/(2*sigma));

% Fokker-Planck evolution generator 
F=v2fplanck(spin_system,parameters); F=inflate(F);

% System trajectory
traj=evolution(spin_system,F,[],A(:),5e-5,200,'trajectory'); traj=full(traj);

% Plotting
kfigure();
for n=1:size(traj,2)
    volplot(reshape(traj(:,n),parameters.npts),[-parameters.dims(1)/2 parameters.dims(1)/2 ...
                                                -parameters.dims(2)/2 parameters.dims(2)/2 ...
                                                -parameters.dims(3)/2 parameters.dims(3)/2]);
    kgrid; drawnow;
end

end

