% Exemplary setup of a 19F Mims ENDOR simulation at 34 GHz using
% Spinach constructor, ENDOR context, and pulse-sequence function style.
%
% March 2024 A. Kehl  (akehl@gwdg.de)
% Spinach architecture migration May 2026 Talos

function endor_kehl_mims_powder()

% Define constants
constants=kehl_constants();

% Spinach spin-system constructor input
T=0.066;
D=6.5;

sys.magnet=1.20521;
sys.isotopes={'E','19F','1H','1H','14N'};
sys.output='hush';
sys.disable={'hygiene'};

inter.zeeman.matrix=cell(1,5);
inter.zeeman.matrix{1}=diag([2.00886,2.00610,2.00211]);
inter.zeeman.matrix{2}=zeros(3,3);
inter.zeeman.eigs=cell(1,5);
inter.zeeman.euler=cell(1,5);
inter.zeeman.eigs{2}=[205,310,358];
inter.zeeman.euler{2}=[87,-14,-29]*pi/180;
inter.zeeman.matrix{3}=zeros(3,3);
inter.zeeman.matrix{4}=zeros(3,3);
inter.zeeman.matrix{5}=zeros(3,3);

inter.coupling.eigs=cell(5,5);
inter.coupling.euler=cell(5,5);
inter.coupling.eigs{1,2}=[2*T,-T,-T]*1e6;
inter.coupling.euler{1,2}=[0,-1,-161]*pi/180;
inter.coupling.eigs{1,5}=[15,11,95.8]*1e6;
inter.coupling.euler{1,5}=[0,0,0];
inter.coupling.eigs{5,5}=diag(remtrace(diag([1.2,0.54,-1.7]*1e6)))';
inter.coupling.euler{5,5}=[0,0,0];
inter.coupling.eigs{2,3}=[2*D,-D,-D]*1e3;
inter.coupling.euler{2,3}=[0,-143,16.9]*pi/180;
inter.coupling.eigs{2,4}=[2*D,-D,-D]*1e3;
inter.coupling.euler{2,4}=[0,-34.3,47.9]*pi/180;

bas.formalism='zeeman-hilb';
bas.approximation='none';

spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
expt=kehl_exp_create(33.7727,12052.1,0.5,0.001,0.3,...
                      [12,2000,12,1000,200000,1000,12]);
expt=kehl_mims_fields(constants,expt,true,[0,0]);
expt=kehl_exp_freq(expt,32,3,0.001);
expt=kehl_exp_angle(expt,163);
expt=kehl_exp_pulse(expt,...
    fullfile(fileparts(mfilename('fullpath')),'kehl_mlr09_12ns_pulse.txt'),false);

% ENDOR context metadata and simulation parameters
parameters.inter=inter;
parameters.freqDomain=true;
parameters.powder=true;
parameters.Nang=50;
parameters.Relax=false;
parameters.Bterm=false;
parameters.endor_spins=[2,3,4];
parameters.epr_spins=5;
parameters.epr_quadrupole_matrix=diag([1.2,0.54,-1.7]*1e6);
parameters.n_spin_systems=1;
parameters.dipolar_pairs=[2,3;2,4];
parameters.expt=expt;

% Actual ENDOR calculation through Spinach-style context and experiment
[endor_amp,endor_amp_conv,x_coords,v_L]=endor_kehl_context(spin_system,@endor_kehl_mims,parameters,'labframe'); %#ok<ASGLU>

% Plotting
sim=endor_amp(:)-endor_amp(1);
sim=sim/min(sim);
plot((x_coords(:))*1e-6,sim);

end

