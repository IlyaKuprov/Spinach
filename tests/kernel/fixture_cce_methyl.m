% Numerical regression fixture for the original disjoint proton model.
% Retains the independent dense-reference comparison in test_cce_examples.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=fixture_cce_methyl()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','1H','1H','1H'};
sys.output='hush';
sys.parallel={'local',1};
sys.tols.cce_atol=1e-8;
sys.tols.cce_rtol=1e-8;
sys.tols.cce_mem=2^28;

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

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:4;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=3;
parameters.cce.approximation='projected';
parameters.cce.expansion='disjoint';
parameters.cce.groups={[2 3 4]};

% Execute the named Hahn echo through the cluster-expansion context
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end


