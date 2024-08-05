% Pulse-acquire experiment during the Diels-Alder cycloaddition 
% of acetylene to butadiene, demonstrating the non-linear kine-
% tics module.
%
% Calculation time: hours, faster on a Tesla A100 GPU.
%
% i.kuprov@soton.ac.uk
% a.acharya@soton.ac.uk

function diels_alder_spec()

% DFT import options
options.min_j=0.5;         % Minimum significant J-coupling, Hz

% Load acetylene      (substance A)
props_a=gparse('acetylene.out');
[sys_a,inter_a]=g2spinach(props_a,{{'H','1H'}},31.8,options);

% Load butadiene      (substance B)
props_b=gparse('butadiene.out');
[sys_b,inter_b]=g2spinach(props_b,{{'H','1H'}},31.8,options);

% Load cyclohexadiene (substance C)
props_c=gparse('cyclohexadiene.out');
[sys_c,inter_c]=g2spinach(props_c,{{'H','1H'}},31.8,options);

% Merge the spin systems
[sys,inter]=merge_inp({sys_a,sys_b,sys_c},{inter_a,inter_b,inter_c});

% Add natural abundance ethanol   (substance D)
sys_d.isotopes={'1H','1H','1H','1H','1H','1H'};
inter_d.zeeman.matrix={1.26, 1.26, 1.26, ...
                       3.69, 3.69, 2.61}*eye(3);
inter_d.coordinates={[]; []; []; []; []; []};
inter_d.coupling.scalar=zeros(6,6);
inter_d.coupling.scalar(1,[4 5])=7.0;
inter_d.coupling.scalar(2,[4 5])=7.0;
inter_d.coupling.scalar(3,[4 5])=7.0;
inter_d.coupling.scalar=num2cell(inter_d.coupling.scalar);
[sys,inter]=merge_inp({sys,sys_d},{inter,inter_d});

% Magnet field
sys.magnet=14.1;

% Initial concentrations, mol/L
A0=1e-2; B0=2e-2; C0=0; D0=0.1;

% Chemical parts and concentrations
inter.chem.parts={1:2, 3:8, 9:16, 17:22};
inter.chem.concs=[A0 B0 C0 D0];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={ 5e-12 ... % Acetylene
             20e-12 ... % Butadiene
             50e-12 ... % Cyclohexadiene
             10e-12};   % Ethanol

% This needs a GPU
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% 2nd order rate constant
k=250.0;  % mol/(L*s)

% 2nd order reaction generator
K=@(t,x)(1i*[-k*x(2)   0        0        0; 
              0       -k*x(1)   0        0;
              0        k*x(1)   0        0;
              0        0        0        0]);

% Time grid (one second)
nsteps=1000; tmax=1.0; dt=tmax/nsteps;
time_axis=linspace(0,tmax,nsteps+1); 
 
% Preallocate trajectory 
x=zeros(4,nsteps+1);

% Define initial concentrations
x(:,1)=[A0 B0 C0 D0]'; 
 
% Run Lie group solver
for n=1:nsteps 
    x(:,n+1)=iserstep(spin_system,K,x(:,n),n*dt,dt,'LG4'); 
end

% Interpolate concentrations as functions of time
A=griddedInterpolant(time_axis,x(1,:),'makima','none');
B=griddedInterpolant(time_axis,x(2,:),'makima','none');
C=griddedInterpolant(time_axis,x(3,:),'makima','none');
D=griddedInterpolant(time_axis,x(4,:),'makima','none');

% Nonlinear kinetics generator
reaction.reactants=[1 2];  % which substances are reactants
reaction.products=3;       % which substances are products
reaction.matching=[1  9;  
                   2 12;   % which spin on the left hand
                   3 15;   % side of the reaction arrow
                   4 16;   % goes into which spin on the
                   5 10;   % right hand side
                   6 11;
                   7 14;
                   8 13];
[GD,GF]=react_gen(spin_system,reaction); 
G=gpuArray(GD{1}+GD{2}+GF{1});

% Set up a pulse-acquire experiment
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.offset=2400;
parameters.sweep=5000;
parameters.npoints=4096;
parameters.zerofill=16384;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Set up time discretisation
dt=1/parameters.sweep;

% Get evolution generators
spin_system=assume(spin_system,'nmr');
H=hamiltonian(spin_system); 
H=frqoffset(spin_system,H,parameters); H=gpuArray(H);
R=relaxation(spin_system); R=gpuArray(R);

% Preallocate the trajectory and get it started
traj=gpuArray.zeros([numel(parameters.rho0) ...
                     parameters.npoints]); 
traj(:,1)=parameters.rho0;

% Run the evolution loop
for n=1:(parameters.npoints-1)

    % Get left and right edge reaction rates
    rr_L=k*A((n-1)*dt)*B((n-1)*dt); rr_R=k*A(n*dt)*B(n*dt);

    % Make left and right edge evolution generators
    F_L=H+1i*R+1i*rr_L*G; F_R=H+1i*R+1i*rr_R*G;

    % Take the time step using Anu's two-point call
    traj(:,n+1)=step(spin_system,{F_L,F_R},traj(:,n),dt);
    
    % Update the user
    report(spin_system,['Non-linear kinetics evolution step ' ...
                        num2str(n) '/' num2str(parameters.npoints-1)]);

end

% Get back to CPU
traj=gather(traj);

% Preallocate the fid
fid=zeros(parameters.npoints,1);

% Run the detection loop
parfor n=1:parameters.npoints
    
    % Localise spin system object
    local_sso=spin_system;

    % Update concentrations
    local_sso.chem.concs=[A((n-1)*dt) ...
                          B((n-1)*dt) ...
                          C((n-1)*dt) ...
                          D((n-1)*dt)]; %#ok<PFBNS>

    % Get chemically weighted detection state
    coil=state(local_sso,'L+','1H','chem');

    % Take the inner product
    fid(n)=coil'*traj(:,n);

end

% Apodization
fid=apodization(fid,'gaussian-1d',10);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);
figure(); plot_1d(spin_system,real(spectrum),parameters); xlim([5.5 6.9]);
figure(); plot_1d(spin_system,real(spectrum),parameters); xlim([2.9 3.0]);

end

