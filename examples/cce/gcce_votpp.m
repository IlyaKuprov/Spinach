% An electron/51V observed multilevel spin coupled to a small proton bath.
% doi:10.1063/5.0278659; illustrative observed multilevel spin, not VOTPP tensors.
% This small illustrative model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=gcce_votpp()
%
% Returns one row of column-vectorised observed density per time point
% in answer, and the cluster-expansion diagnostics structure in diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=gcce_votpp()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0;
sys.isotopes={'E','51V','1H'};
sys.output='hush';
sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,2}=diag([15 15 50])*1e6;
inter.coupling.matrix{1,3}=diag([2000 3000 7000]);
inter.coupling.matrix{2,3}=diag([100 100 -200]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*1e5;

% Basis formalism for retained clusters
bas.formalism='zeeman-hilb';

% Experiment parameters
parameters.observed=[1 2];
parameters.transition=struct('spins',[1 2],'ms',[-0.5 0.5;3.5 3.5]);
parameters.output='density';
parameters.times=linspace(0,200e-6,31);
parameters.cce.approximation='generalized';
parameters.cce.order=1;
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


