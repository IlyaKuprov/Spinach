% Optomechanical sideband transfer of a phonon Fock state into a
% driven cavity. A red-detuned coherent drive on the cavity acti-
% vates the beam-splitter part of the radiation pressure coupling,
% and the second Fock state of the mechanical oscillator is cohe-
% rently exchanged with the cavity field. Model and parameters
% from the propagation test set of QuantumPropagators.jl; all
% quantities are in the dimensionless units of the source, with
% input values divided by 2*pi to cancel the Hz convention.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function optomech_sideband()

% Magnet field
sys.magnet=0;

% Cavity with five and phonon mode with eleven Fock levels
sys.isotopes={'C5','V11'};

% Mode frequencies, cavity in the red-detuned drive rotating frame
inter.modes.frqs={10/(2*pi) 10/(2*pi)};

% Radiation pressure coupling, cavity number times phonon coordinate
inter.modes.longitudinal=cell(2,2);
inter.modes.longitudinal{1,2}=-sqrt(2)/(2*pi);

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'labframe'));

% Coherent drive on the cavity
H=H+2*(operator(spin_system,'C',1)+operator(spin_system,'A',1));

% Number operators of the two modes
Nc=operator(spin_system,'N',1);
Nm=operator(spin_system,'N',2);

% Empty cavity and second phonon Fock state
rho=state(spin_system,{'BL1','BL3'},{1,2});

% Time grid of the source calculation
nsteps=250; dt=0.2; time_axis=dt*(0:nsteps);

% Propagator over one time step
P=propagator(spin_system,H,dt);

% Preallocate occupation trajectories
n_cav=zeros(nsteps+1,1); n_mech=zeros(nsteps+1,1);

% Propagate and record mode occupations
for n=1:(nsteps+1)
    n_cav(n)=real(trace(Nc*rho));
    n_mech(n)=real(trace(Nm*rho));
    rho=P*rho*P';
end

% Validate the trace preservation
if abs(trace(rho)-1)>1e-6
    error('propagation is not trace-preserving.');
end

% Validate the initial phonon occupation
if abs(n_mech(1)-2)>1e-12
    error('initial phonon occupation is not two.');
end

% Validate the sideband transfer into the cavity
if max(n_cav)<1.5
    error('sideband transfer into the cavity did not occur.');
end

% Plot the mode occupation dynamics
kfigure(); plot(time_axis,[n_cav n_mech],'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, dimensionless');
kylabel('mode occupation numbers');
ktitle('optomechanical sideband transfer');
klegend({'cavity','mechanical mode'},'Location','Best');

end

