% Non-linear reaction kinetics in a situation when there is
% no hydrodynamics, diffusion, or spin dynamics. This is in-
% tended as a stepping stone to the more complicated cases
% in the same directory of the Spinach example set. Explicit reaction
% records act on five spin-free concentration coordinates. Additive
% product arrival shares its unit source equally between reactants.
%
% Calculation time: seconds.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function plain_reaction()

% Trace a ghost seed to five concentration pools, independent of solvent spins
sys.magnet=0; sys.isotopes={'G'};
inter.chem.parts={1,[],[],[],[]};
inter.chem.concs=[0.6 0.5 0.0 0.0 18.1];

% Competing endo and exo cycloadditions, rates in L/(mol*s)
inter.chem.reactions={...
    struct('reactants',[1 2],'products',3,'matching',zeros(0,2),'rate',0.5),...
    struct('reactants',[1 2],'products',4,'matching',zeros(0,2),'rate',0.1)};
bas.formalism='sphten-liouv';
bas.approximation={'none','none','none','none','none'};

% Build the concentration-only system and compile its reaction network
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=kill_spin(spin_system,1);
K=kinetics(spin_system);

% Kinetic time grid, 20 seconds
nsteps=200; tmax=20; dt=tmax/nsteps;
time_axis=linspace(0,tmax,nsteps+1); 
 
% Preallocate concentration trajectory
x=zeros(5,nsteps+1);

% Initial concentrations, mol/L
x(:,1)=unit_state(spin_system);

% Concentration dynamics
for n=1:nsteps 
    x(:,n+1)=step(spin_system,{@(t,y)1i*K(t,y),n*dt,'LG4'},x(:,n),dt);
end

% Plot concentrations, excluding solvent
kfigure(); plot(time_axis,real(x(1:4,:))); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'cyclopentadiene','acrylonitrile', ...
         'endo-norbornene carbonitrile',...
         'exo-norbornene carbonitrile'}, ...
         'Location','Best'); drawnow;

end

