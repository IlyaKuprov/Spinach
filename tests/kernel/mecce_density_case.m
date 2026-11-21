% An observed spin coupled to a phenomenologically relaxing spin-half bath.
% The independent dense Liouvillian comparison is in test_cce_examples.
% The supplied rate is not a prediction of microscopic T1. See
% https://doi.org/10.1103/PhysRevLett.132.250401. Syntax:
%
%            [answer,diagnostics]=mecce_density_case()
%
% Returns one row of column-vectorised observed density per time point
% in answer, and the cluster-expansion diagnostics structure in diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=mecce_density_case()

% Standard Spinach system and anisotropic hyperfine tensor in Hz
sys.magnet=0;
sys.isotopes={'E','13C'};
sys.output='hush';
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([1100 900 2200]);

% Basis formalism for retained clusters
bas.formalism='zeeman-liouv';

% Experiment parameters
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[0.5 -0.5]);
parameters.output='density';
parameters.times=linspace(0,200e-6,31);
parameters.bath_jumps=struct('support',2,'operator','L+','rate',3000);
parameters.cce.approximation='generalized';
parameters.cce.order=1;
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);

% Run the named experiment in the cluster-expansion context
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@free_evolution,parameters,'labframe');

end


