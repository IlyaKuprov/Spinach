% Effective electron with anisotropic 183W couplings: reduced rare-earth echo.
% doi:10.1126/sciadv.abj9786; CIF, g-tensor, and exclusion policy are gated.
% The deliberately small model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=cce_erbium_cawo4_small()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_erbium_cawo4_small()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','183W','183W'};
sys.output='hush';
sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[8 0 2];[0 10 -1]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;3000 1200 2000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;-900 400 1000].';
inter.coupling.matrix{2,3}=diag([80 80 -160]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Basis policy for retained clusters
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:3;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=2;
parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);

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


