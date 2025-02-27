% Non-linear reaction kinetics in a situation when there is
% no hydrodynamics, diffusion, spin dynamics, or relaxation, 
% but magnetic observables are still being detected. This is
% intended as a stepping stone to the more complicated cases
% in the same directory of the Spinach example set.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function plain_kinetics()

% Import Diels-Alder cycloaddition
[sys,inter,bas]=dac_reaction();

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;

% This needs a GPU
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rate constants, mol/(L*s)
k1=250; % towards exo  
k2=50;  % towards endo

% Cycloaddition reaction generator, including solvent
K=@(t,x)(1i*[-k1*x(2)-k2*x(2)  0                0      0     0;      
              0               -k1*x(1)-k2*x(1)  0      0     0;           
              0                k1*x(1)          0      0     0;
              0                k2*x(1)          0      0     0; 
              0                0                0      0     0]);        
                     
% Time grid (one second)
nsteps=200; tmax=0.2; dt=tmax/nsteps;
time_axis=linspace(0,tmax,nsteps+1); 
 
% Preallocate concentration trajectory
x=zeros(5,nsteps+1);

% Initial concentrations, mol/L
x(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]; 

% Stage 1: concentration dynamics
for n=1:nsteps 
    x(:,n+1)=iserstep(spin_system,K,x(:,n),n*dt,dt,'LG4'); 
end

% Plot concentrations, excluding solvent
figure(); scale_figure([1.81 0.85]);
subplot(1,2,1); plot(time_axis,real(x(1:4,:))); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'cyclopentadiene','acrylonitrile', ...
         'endo-norbornene carbonitrile',...
         'exo-norbornene carbonitrile'}, ...
         'Location','northeast'); drawnow;

% Interpolate concentrations as functions of time
A=griddedInterpolant(time_axis,x(1,:),'makima','none');
B=griddedInterpolant(time_axis,x(2,:),'makima','none');
C=griddedInterpolant(time_axis,x(3,:),'makima','none');
D=griddedInterpolant(time_axis,x(4,:),'makima','none');
E=griddedInterpolant(time_axis,x(5,:),'makima','none');

% Nonlinear kinetics generator
reaction{1}.reactants=[1 2];  % cyclopentadiene and acrylonitrile
reaction{1}.products=3;       % into endo-norbornene carbonitrile
reaction{1}.matching=[1 12; 2 17; 3 18; 4 16; 5 10; 6 11; 7 14; 8 15; 9 13];
reaction{2}.reactants=[1 2];  % cyclopentadiene and acrylonitrile
reaction{2}.products=4;       % into endo-norbornene carbonitrile
reaction{2}.matching=[1 21; 2 26; 3 27; 4 25; 5 19; 6 20; 7 23; 8 24; 9 22]; 

% Reaction generators
G1=react_gen(spin_system,reaction{1});
G2=react_gen(spin_system,reaction{2});

% Get concentration-weighted initial condition
eta=A(0)*state(spin_system,'L+',spin_system.chem.parts{1})+...
    B(0)*state(spin_system,'L+',spin_system.chem.parts{2})+...
    C(0)*state(spin_system,'L+',spin_system.chem.parts{3})+...
    D(0)*state(spin_system,'L+',spin_system.chem.parts{4})+...
    E(0)*state(spin_system,'L+',spin_system.chem.parts{5});

% Preallocate the trajectory 
traj=zeros([numel(eta) nsteps+1]); traj(:,1)=eta;

% No spin evolution here
H=sparse(0); R=sparse(0);

% Stage 2: nuclear spin dynamics
for n=1:nsteps

    % Build the left interval edge composite evolution generator
    F_L=H+1i*R+1i*k1*G1{1}*B(time_axis(n))...   % Reaction 1 from substance A
              +1i*k1*A(time_axis(n))*G1{2}...   % Reaction 1 from substance B
              +1i*k2*G2{1}*B(time_axis(n))...   % Reaction 2 from substance A
              +1i*k2*A(time_axis(n))*G2{2};     % Reaction 2 from substance B

    % Build the right interval edge composite evolution generator
    F_R=H+1i*R+1i*k1*G1{1}*B(time_axis(n+1))... % Reaction 1 from substance A
              +1i*k1*A(time_axis(n+1))*G1{2}... % Reaction 1 from substance B
              +1i*k2*G2{1}*B(time_axis(n+1))... % Reaction 2 from substance A
              +1i*k2*A(time_axis(n+1))*G2{2};   % Reaction 2 from substance B

    % Take the time step using the two-point Lie quadrature
    traj(:,n+1)=step(spin_system,{F_L,F_R},traj(:,n),time_axis(n+1)-time_axis(n));

end

% Detect particular spins 
coils=[state(spin_system,{'L+'},{1})  ...
       state(spin_system,{'L+'},{7})  ...
       state(spin_system,{'L+'},{10}) ...
       state(spin_system,{'L+'},{19})];
x=real((coils'*coils)\(coils'*traj));

% Plot observables, excluding solvent
subplot(1,2,2); plot(time_axis,x); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('observable, mol/L');
klegend({'cyclopentadiene','acrylonitrile', ...
         'endo-norbornene carbonitrile',...
         'exo-norbornene carbonitrile'}, ...
         'Location','northeast'); drawnow;

end

