% Exemplary setup of a 19F time-domain ENDOR simulation at 34 GHz using
% Spinach constructor, ENDOR context, and pulse-sequence function style.
%
% March 2024 A. Kehl  (akehl@gwdg.de)
% Spinach architecture migration May 2026 Talos

function endor_kehl_time_powder()

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
inter.zeeman.matrix{3}=zeros(3,3);
inter.zeeman.matrix{4}=zeros(3,3);
inter.zeeman.matrix{5}=zeros(3,3);

inter.coupling.matrix=cell(5,5);
inter.coupling.matrix{1,2}=local_tensor([2*T,-T,-T]*1e6,[161,1,0]);
inter.coupling.matrix{1,5}=local_tensor([15,11,95.8]*1e6,[0,0,0]);
inter.coupling.matrix{5,5}=local_traceless_tensor([1.2,0.54,-1.7]*1e6,[0,0,0]);
inter.coupling.matrix{2,3}=local_tensor([2*D,-D,-D]*1e3,[-16.9,143.0,0]);
inter.coupling.matrix{2,4}=local_tensor([2*D,-D,-D]*1e3,[-47.9,34.3,0]);

bas.formalism='zeeman-hilb';
bas.approximation='none';

spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
expt=kehl_exp_create(33.7727,12052.1,0.5,0.24,144.24,...
                      [12,1600,12,1000,6000,1000,6000,344240,12]);
expt=kehl_time_fields(constants,expt,true,[0,0]);
expt=kehl_exp_freq(expt,32,3,0.001);
expt=kehl_exp_pulse(expt,...
    fullfile(fileparts(mfilename('fullpath')),'kehl_mlr09_12ns_pulse.txt'),false);

% ENDOR context metadata and simulation parameters
parameters.inter=inter;
parameters.freqDomain=true;
parameters.powder=true;
parameters.Nang=50;
parameters.endor_spins=[2,3,4];
parameters.epr_spins=5;
parameters.epr_quadrupole_matrix=local_tensor([1.2,0.54,-1.7]*1e6,[0,0,0]);
parameters.n_spin_systems=1;
parameters.dipolar_pairs=[2,3;2,4];
parameters.expt=expt;
parameters.time_domain=true;

% Actual ENDOR calculation through Spinach-style context and experiment
[endor_amp,endor_amp_conv,x_coords,v_L]=endor_kehl_context(spin_system,@endor_kehl_time,parameters,'labframe'); %#ok<ASGLU>

% Plotting
sim=endor_amp(:)-endor_amp(1);
sim=sim/min(sim);
plot((x_coords(:))*1e-6,sim);

end

function T=local_traceless_tensor(values,eulers_deg)
values=values-mean(values);
T=local_tensor(values,eulers_deg);
end

function T=local_tensor(values,eulers_deg)
alpha=eulers_deg(1)*pi/180;
beta=eulers_deg(2)*pi/180;
gamm=eulers_deg(3)*pi/180;
R=zeros(3,3);
R(1,1)=cos(beta)*cos(alpha)*cos(gamm)-sin(alpha)*sin(gamm);
R(1,2)=cos(beta)*sin(alpha)*cos(gamm)+cos(alpha)*sin(gamm);
R(1,3)=-sin(beta)*cos(gamm);
R(2,1)=-cos(beta)*cos(alpha)*sin(gamm)-sin(alpha)*cos(gamm);
R(2,2)=-cos(beta)*sin(alpha)*sin(gamm)+cos(alpha)*cos(gamm);
R(2,3)=sin(beta)*sin(gamm);
R(3,1)=sin(beta)*cos(alpha);
R(3,2)=sin(beta)*sin(alpha);
R(3,3)=cos(beta);
T=R*diag(values)*R';
end

