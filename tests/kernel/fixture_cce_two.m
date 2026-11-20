% Finite four-spin gCCE density fixture for the full-Hilbert reference.
% Returns the normalised vectorised observed densities and CCE diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=fixture_cce_two()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0;
sys.isotopes={'E','E','13C','13C'};
sys.output='hush';
sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,2}=diag([700 500 1200]);
inter.coupling.matrix{1,3}=diag([1700 900 2000]);
inter.coupling.matrix{2,4}=diag([1300 800 1600]);
inter.coupling.matrix{3,4}=diag([100 100 -200]);

% Basis formalism for retained clusters
bas.formalism='zeeman-hilb';

% Experiment parameters
parameters.observed=[1 2];
parameters.transition=struct('spins',[1 2],'ms',[0.5 -0.5;0.5 -0.5]);
parameters.output='density';
parameters.times=linspace(0,200e-6,31);
parameters.cce.approximation='generalized';
parameters.cce.order=2;
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);

% Run the named experiment in the cluster-expansion context
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@free_evolution,parameters,'labframe');

end


