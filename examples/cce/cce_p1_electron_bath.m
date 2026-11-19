% Sparse electron bath with an explicit internal nitrogen hyperfine group.
% doi:10.1103/PhysRevB.110.134106; one P1-like group, not a concentration fit.
% The deliberately small model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=cce_p1_electron_bath()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_p1_electron_bath()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','E','14N'};
sys.output='hush';
sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[35 0 0];[35 0 1]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=diag([2000 2000 -4000]);
inter.coupling.matrix{2,3}=diag([81 81 114])*1e6;
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;400 0 200].';

% Declare observed-Sz projection for this projected electron-bath model
for n=2:numel(sys.isotopes)

    % Remove transverse observed-spin coupling components
    if ~isempty(inter.coupling.matrix{1,n})
        inter.coupling.matrix{1,n}(:,1:2)=0;
    end

end

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Local Hilbert basis for retained clusters
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:3;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=1;
parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-3,'rtol',1e-3,...
                               'max_memory_bytes',2^28);

% Keep the electron and its nitrogen in the same physical group
parameters.cce.expansion='partition';
parameters.cce.groups={[2 3]};
parameters.cce.partition=struct('mode','user','max_dimension',6);

% Execute the named Hahn experiment through the cluster context
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Plot the measured complex echo amplitude
kfigure();
plot(1e6*parameters.times,real(answer),'b-',...
     1e6*parameters.times,imag(answer),'r--');
xlabel('total echo time / microseconds');
ylabel('normalised echo amplitude E(2τ)');
legend('real','imaginary');
kgrid();

end


