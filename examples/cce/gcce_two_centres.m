% Two electron centres: density evolution from a selected transition.
% doi:10.1063/5.0303171 motivates the category; geometry is illustrative.
% This small illustrative model is not a published material reproduction.
% Independent density and negativity checks are in test_cce_examples.
% Syntax:
%
%             [answer,diagnostics]=gcce_two_centres()
%
% Returns one row of column-vectorised observed density per time point
% in answer, and the cluster-expansion diagnostics structure in diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=gcce_two_centres()

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

% Plot the physical observed-spin response
populations=real(answer(:,1:sqrt(size(answer,2))+1:end));
kfigure();
plot(1e6*parameters.times,populations);
xlabel('time / microseconds');
ylabel('observed populations');
kgrid();

end


