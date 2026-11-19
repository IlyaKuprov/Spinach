% A three-proton methyl-like group with coherent cyclic exchange.
% doi:10.5194/mr-6-15-2025 motivates grouping, not these chosen tensors.
% The deliberately small model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=cce_methyl_molecular()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_methyl_molecular()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','1H','1H','1H'};
sys.output='hush';
sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[5 1 0];[5 -0.5 0.866];[5 -0.5 -0.866]};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;7000 0 5000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;5000 0 4000].';
inter.coupling.matrix{1,4}=[0 0 0;0 0 0;6000 0 3000].';
inter.coupling.matrix{2,3}=2000*eye(3);
inter.coupling.matrix{2,4}=2000*eye(3);
inter.coupling.matrix{3,4}=2000*eye(3);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Basis policy for retained clusters
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:4;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=3;
parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);
parameters.cce.expansion='disjoint';
parameters.cce.groups={[2 3 4]};

% Execute the named Hahn echo through the cluster-expansion context
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Plot the complex echo amplitude
figure();
plot(1e6*parameters.times,real(answer),'b-',...
     1e6*parameters.times,imag(answer),'r--');
xlabel('total echo time / microseconds');
ylabel('normalised echo amplitude E(2τ)');
legend('real','imaginary');
kgrid();

end


