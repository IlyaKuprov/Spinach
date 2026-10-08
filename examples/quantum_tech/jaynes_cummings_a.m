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

% Cavity resonant with the electron
e_frq=-sys.magnet*spin('E')/(2*pi);
inter.modes.frqs={[] e_frq};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=2.828e6;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotating frame Hamiltonian
spin_system=assume(spin_system,'cavity');
H_JC=hamiltonian(spin_system);

% Electron detuning operator
Ez=operator(spin_system,{'Lz'},{1});

% Locate the one-excitation manifold
spin_exc=state(spin_system,{'ZL2','BL1'},{1,2});
cav_exc=state(spin_system,{'ZL1','BL2'},{1,2});
one_quant=speye(size(H_JC,1));
one_quant=one_quant(:,[find(diag(spin_exc)>0.5)...
                       find(diag(cav_exc)>0.5)]);

% Detuning range
delta=2*pi*linspace(-15e6,15e6,100);

% Eigenvalue array
eig_array=zeros(2,100);

% Loop over detunings
for n=1:numel(delta)

    % Make the Hamiltonian
    H=delta(n)*Ez+H_JC; H=(H+H')/2;

    % Diagonalise the one-photon manifold
    eig_array(:,n)=eig(full(one_quant'*H*one_quant));

end

% Plot one-photon case
kfigure(); plot(1e-6*delta/(2*pi),...
                1e-6*eig_array/(2*pi));
axis tight; kxlabel('detuning, MHz');
kylabel('energy levels, MHz'); kgrid;

end

