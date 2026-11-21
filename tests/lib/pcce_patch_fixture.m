% Partitioned CCE for a strongly coupled electron pair in a planar layer.
% Four-spin finite-patch regression fixture; not a material reproduction.
% Syntax:
%
%             [answer,diagnostics]=pcce_patch_fixture()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=pcce_patch_fixture()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','E','E','E'};
sys.output='hush';
sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[20 0 0];[21 0 0];[0 30 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;-12000 0 6000].';
inter.coupling.matrix{1,4}=[0 0 0;0 0 0;4000 0 3000].';
inter.coupling.matrix{2,3}=diag([2e5 2e5 -4e5]);
inter.coupling.matrix{3,4}=diag([300 300 -600]);

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
parameters.cce.bath_spins=2:4;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=2;
parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-5,'rtol',1e-5,...
                               'max_memory_bytes',2^28);
parameters.cce.expansion='partition';
parameters.cce.groups={[2 3],4};
parameters.cce.partition=struct('mode','user','max_dimension',4);

% Execute the named Hahn experiment through the cluster context
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end


