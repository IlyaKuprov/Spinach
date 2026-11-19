% Spin-3/2 vacancy-like centre with mixed spin-half nuclear environment.
% doi:10.1103/PhysRevA.109.022603; lattice and transition geometry gate.
% This small illustrative model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=gcce_highspin_siv()
%
% Returns one row of column-vectorised observed density per time point
% in answer, and the cluster-expansion diagnostics structure in diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=gcce_highspin_siv()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0;
sys.isotopes={'E4','13C','29Si'};
sys.output='hush';
sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,1}=diag([-1 -1 2])*1e4;
inter.coupling.matrix{1,2}=diag([1700 1300 2400]);
inter.coupling.matrix{1,3}=diag([800 900 1500]);
inter.coupling.matrix{2,3}=diag([100 100 -200]);

% Basis formalism for retained clusters
bas.formalism='zeeman-hilb';

% Experiment parameters
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[-0.5 0.5]);
parameters.output='density';
parameters.times=linspace(0,200e-6,31);
parameters.cce.approximation='generalized';
parameters.cce.order=2;
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);

% Spin system and free evolution
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@free_evolution,parameters,'labframe');

% Plot populations in the observed tensor-product order
populations=real(answer(:,1:sqrt(size(answer,2))+1:end));
kfigure();
plot(1e6*parameters.times,populations);
xlabel('time / microseconds');
ylabel('observed populations');
kgrid();

end


