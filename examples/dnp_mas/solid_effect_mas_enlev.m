% A MAS DNP simulation performed as described in Fred Mentink-
% Vigier's paper (Spinach rotation conventions are different):
%
%         http://dx.doi.org/10.1016/j.jmr.2015.07.001
%
% Energy level diagram as a function of the rotor phase.
%
% Calculation time: milliseconds
%
% ilya.kuprov@weizmann.ac.il

function solid_effect_mas_enlev()

% Magnet field
sys.magnet=9.403;

% Spin specification
sys.isotopes={'E','1H'};

% Interactions
inter.zeeman.eigs{1}=[2.00614  2.00194 2.00988];
inter.zeeman.euler{1}=pi*[253.6 105.1 123.8]/180;
inter.zeeman.eigs{2}=[0.00 0.00 0.00];
inter.zeeman.euler{2}=[0.00 0.00 0.00];
inter.coordinates={[0.00 0.00 0.00],...
                   [0.00 0.00 3.00]};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Stack generation parameters
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rframes={};
parameters.orientation=[0 0 0];
parameters.spins={'E','1H'};
parameters.masframe='magnet';
parameters.offset=[0 0];
parameters.max_rank=200;

% Stack generation
H=rotor_stack(spin_system,parameters,'esr');

% Stack diagonalization
energies=zeros(size(H{1},1),numel(H));
parfor n=1:numel(H)
    energies(:,n)=sort(real(eig(H{n})));
end

% Plotting
phi_axis=linspace(0,2*pi,numel(H)+1)';
phi_axis(end)=[]; kfigure(); 
plot(phi_axis,energies','b-'); kgrid;
axis tight; kxlabel('rotor phase, rad'); 
kylabel('level energy, rad/s');

end

