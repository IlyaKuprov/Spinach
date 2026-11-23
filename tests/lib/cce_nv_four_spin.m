% Projected NV Hahn echo for the pinned PyCCE four-nucleus numerical fixture.
% Source: MICCoMpy/PyCCE d74d3cb683bf3caa07231654314e554f440833b4,
% https://doi.org/10.1002/adts.202100254. Explicit numerical fixture arrays
% are public-interface inputs, not copied implementation source. Syntax:
%
%                         [answer,diagnostics]=cce_nv_four_spin()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_nv_four_spin()

% Load fixed physical inputs; choose the execution resources explicitly
[sys,inter,fixture]=cce_fixture('nv_four_c13','6f7282b6');
sys.output='hush';
sys.parallel={'local',1};
sys.tols.cce_atol=1e-3;
sys.tols.cce_rtol=1e-3;
sys.tols.cce_mem=2^28;

% Declare the reference's observed-Sz projected approximation before propagation
for n=2:5
    inter.coupling.matrix{1,n}(:,1:2)=0;
end

% Spinach housekeeping
spin_system=create(sys,inter);

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:5;
parameters.transition.spins=1;
parameters.transition.ms=[-1 1];
parameters.times=fixture.times;
parameters.cce.order=2;
parameters.cce.approximation='projected';

% Execute the named Hahn experiment through the cluster context
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end


