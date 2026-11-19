% Low-field clock-like spin-one centre with transverse ZFS and a nuclear bath.
% doi:10.1103/PRXQuantum.2.010311; source D/E and atomistic bath gate.
% This small illustrative model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=gcce_sic_clock()
%
% Returns one row of column-vectorised observed density per time point
% in answer, and the cluster-expansion diagnostics structure in diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=gcce_sic_clock()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0;
sys.isotopes={'E3','13C','29Si'};
sys.output='hush';
sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,1}=diag([-1.0e4+2e3 -1.0e4-2e3 2e4]);
inter.coupling.matrix{1,2}=diag([2200 1800 3000]);
inter.coupling.matrix{1,3}=diag([900 1200 1700]);
inter.coupling.matrix{2,3}=diag([120 120 -240]);

% Basis formalism for retained clusters
bas.formalism='zeeman-hilb';

% Experiment parameters
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[-1 1]);
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


