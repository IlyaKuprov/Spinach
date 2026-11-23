% Planar high-spin 11B/14N bath with explicit quadrupolar tensors.
% For material-specific Hahn echoes, see Haykal et al., Fig. 4a,b,
% doi:10.1038/s41467-022-31743-0. This is a two-nucleus model, not a
% reproduction of the bulk gCCE3 calculation.
% Syntax:
%
%             [answer,diagnostics]=cce_hbn_small()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_hbn_small()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','11B','14N'};
sys.output='hush';
sys.parallel={'local',1};
sys.tols.cce_atol=1e-8;
sys.tols.cce_rtol=1e-8;
sys.tols.cce_mem=2^28;

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[1.45 0 0];[-0.725 1.2557 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;40000 0 20000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;20000 0 -10000].';
inter.coupling.matrix{2,3}=diag([500 500 -1000]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
inter.coupling.matrix{3,3}=diag([-1 -1 2])*5e4;

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:3;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=2;
parameters.cce.approximation='projected';

% Execute the named Hahn echo through the cluster-expansion context
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end


