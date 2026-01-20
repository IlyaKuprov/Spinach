% Jaynes-Cummings coupling between a spin and an electromagnetic
% cavity mode with five population numbers included. The avoided
% crossing in the one-photon energy level splitting of the mode
% is plotted with respect to the detuning.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function jaynes_cummings_a()

% Magnet field
sys.magnet=0.33;

% System
sys.isotopes={'E','C5'};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Larmor and cavity energy operators
Ez=operator(spin_system,{'Lz'},{1});
N=operator(spin_system,{'N'},{2});
U=unit_oper(spin_system);

% Electron larmor and cavity frequency
omega_c=-sys.magnet*spin('E');

% Jaynes-Cummings term
g=2*pi*2.828e6;
H_JC=g*(operator(spin_system,{'L+','A'},{1,2})+...
        operator(spin_system,{'L-','C'},{1,2}));

% Detuning range
delta=2*pi*linspace(-15e6,15e6,100);

% Eigenvalue array
eig_array=zeros(10,100);

% Loop over detunings
for n=1:numel(delta)

    % Make the Hamiltonian
    H=delta(n)*Ez+omega_c*Ez;
    H=H+omega_c*(N+U/2);
    H=H+H_JC; H=(H+H')/2;

    % Diagonalise the Hamiltonian
    eig_array(:,n)=eig(H);

end

% Plot one-photon case
kfigure(); plot(1e-6*delta/(2*pi),...
                1e-6*eig_array(2:3,:)/(2*pi));
axis tight; kxlabel('detuning, MHz'); 
kylabel('energy levels, MHz'); kgrid;

end

