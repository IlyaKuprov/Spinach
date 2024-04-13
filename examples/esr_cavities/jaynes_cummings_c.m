% An exchange-coupled two-electron system with the electrons 
% having independent Jaynes-Cummings couplings to the same
% mode of an electromagnetic cavity. A time-domain simulati-
% on starting with transverse spin magnetisation and empty 
% cavity mode. Detected on the Lx operator of the spin and
% magnetic field operator of the cavity mode.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function jaynes_cummings_c()

% Magnet field
sys.magnet=0.33;

% System
sys.isotopes={'E','E','C5'};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Interaction parameters
g1=2*pi*2.828e6; % J-C, first electron
g2=2*pi*2.728e6; % J-C, second electron
J=2*pi*5e6;      % Exchange coupling

% Exchange coupling Hamiltonian
H_EX=J*(operator(spin_system,{'Lx','Lx'},{1,2})+...
        operator(spin_system,{'Ly','Ly'},{1,2})+...
        operator(spin_system,{'Lz','Lz'},{1,2}));

% Jaynes-Cummings Hamiltonian
H_JC=g1*(operator(spin_system,{'L+','An'},{1,3})+ ...
         operator(spin_system,{'L-','Cr'},{1,3}))+...
     g2*(operator(spin_system,{'L+','An'},{2,3})+ ...
         operator(spin_system,{'L-','Cr'},{2,3}));

% Detuning operators and values
delta1=2*pi*5e6; delta2=2*pi*5e6;
Ez1=operator(spin_system,{'Lz'},{1});
Ez2=operator(spin_system,{'Lz'},{2});

% Rotating frame Hamiltonian
H=delta1*Ez1+delta2*Ez2+H_JC+H_EX;

% Initial state Sx on electrons
% in an empty cavity mode
rho=state(spin_system,{'Lx','Em'},{1,3})+...
    state(spin_system,{'Lx','Em'},{2,3});

% Detection state, spins
S=state(spin_system,{'Lx'},{1})+...
  state(spin_system,{'Lx'},{2});

% Detection state, cavity mode
B=(state(spin_system,{'Cr'},{2})-...
   state(spin_system,{'An'},{2}))/2i;

% Run evolution
traj_s=evolution(spin_system,H,S,rho,10e-9,250,'observable');
traj_c=evolution(spin_system,H,B,rho,10e-9,250,'observable');

% Plot the observables
time_axis=linspace(0,2.5,251); % us
figure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(time_axis,real(traj_s));
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$S_{X}$ of the spin');
ktitle('spin dynamics');
subplot(1,2,2); plot(time_axis,real(traj_c));
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$(a^{+}-a)/2i$ of the mode');
ktitle('cavity mode dynamics');

end

