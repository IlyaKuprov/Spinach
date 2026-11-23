% Numerical regression fixture for the original three-node deuteron model.
% Retains the independent dense-reference comparison in test_cce_examples.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=fixture_cce_powder()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','2H','2H'};
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
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;40000 0 20000];
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;20000 0 -10000];
inter.coupling.matrix{2,3}=diag([500 500 -1000]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
inter.coupling.matrix{3,3}=diag([-1 -1 2])*5e4;

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Experiment parameters
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[0.5 -0.5]);
parameters.preparation='transverse';
parameters.times=linspace(0,200e-6,41);

% Cluster expansion controls
parameters.cce.bath_spins=2:3;
parameters.cce.order=2;
parameters.cce.approximation='generalized';

% Average complete Hahn echo expansions over the quadrature
parameters.grid=struct('alphas',[0 0 pi/2],'betas',[0 pi/2 pi/2],...
                       'gammas',[0 0 0],'weights',ones(3,1)/3);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end


