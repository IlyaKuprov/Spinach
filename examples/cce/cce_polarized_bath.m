% Nuclear polarisation with a sampled product-state bath.
% A two-proton toy bath illustrates thermal population control; it is not
% the molecular geometry in https://doi.org/10.1039/d6cp01591c.
% The exhaustive reference is evaluated in test_cce_examples.
% Syntax:
%
%             [answer,diagnostics]=cce_polarized_bath()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_polarized_bath()

% Standard Spinach magnet, isotope list, and illustrative Hz tensors
sys.magnet=0.01;
sys.isotopes={'E','1H','1H'};
sys.output='hush';
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;9000 0 12000].';
inter.coupling.matrix{2,3}=diag([400 400 -800]);

% Basis formalism for retained clusters
bas.formalism='zeeman-hilb';

% Experiment parameters
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[0.5 -0.5]);
parameters.times=linspace(0,200e-6,21);
parameters.cce.approximation='projected';
parameters.cce.order=2;
parameters.cce.bath_state=struct('mode','sampled','basis','uniform',...
                                 'probabilities',{{[0.8 0.2],[0.8 0.2]}},'samples',64,'seed',8800);
parameters.cce.averaging='after_expansion';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);

% Run the named experiment in the cluster-expansion context
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Plot the physical observed-spin response
kfigure();
plot(1e6*parameters.times,real(answer),'b-');
xlabel('total echo time / microseconds');
ylabel('normalised echo amplitude E(2τ)');
legend('64 sampled states');
kgrid();

end


