% Basic two-transmon system with Duffing model inter-
% actions and a flip-flop coupling.
%
% ilya.kuprov@weizmann.ac.il

function two_transmons()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3','T5'};

% Interactions
inter.duffing.offset={1e6 2e6};
inter.duffing.anharm={1e3 2e3};
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=1e3;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Assumptions
spin_system=assume(spin_system,'duffing');

end

