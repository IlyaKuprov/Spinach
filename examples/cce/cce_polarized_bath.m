% Nuclear polarization and exhaustive versus sampled bath estimators.
% A two-proton toy bath illustrates thermal population control; it is not
% the molecular geometry in https://doi.org/10.1039/d6cp01591c.
% Syntax: result=cce_polarized_bath()
%
% ilya.kuprov@weizmann.ac.il

function result=cce_polarized_bath()

% Standard Spinach magnet, isotope list and illustrative Hz tensors
sys.magnet=0.01; sys.isotopes={'E','1H','1H'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000];
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;9000 0 12000];
inter.coupling.matrix{2,3}=diag([400 400 -800]);
spin_system=create(sys,inter);

% Hahn coherence at full order for the two bath spins
parameters.central_spins=1; parameters.bath_spins=[2 3];
parameters.method='cce'; parameters.expansion='standard';
parameters.quantity='coherence'; parameters.order=2;
parameters.central.alpha=[1;0]; parameters.central.beta=[0;1];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete'); parameters.times=linspace(0,200e-6,21);
parameters.sequence=struct('type','ideal','target',1,'axis','x',...
                           'angle',pi,'placement','fraction','at',0.5);
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Exhaustively average a polarized product density without sampling error
population=[0.8;0.2]; parameters.bath_state=struct('mode','explicit_product',...
    'rho',{{diag(population),diag(population)}});
exact=cce(spin_system,parameters);

% Sample the same distribution with deterministic Threefry state IDs
parameters.bath_state=struct('mode','sampled','basis',{{eye(2),eye(2)}},...
    'probabilities',{{population,population}},'samples',64,'seed',8800);
parameters.averaging='after_expansion'; result=cce(spin_system,parameters);
result.exhaustive=exact.value; result.sampling_error=max(abs(result.value-exact.value));
fprintf('Polarized bath: 64-state sample/exhaustive difference %.3g.\n',...
        result.sampling_error);
fprintf('Published point-centre T2: 9.06,10.30,13.77 microseconds; no toy fit.\n');

% Display finite-sample uncertainty rather than claiming exact agreement
kfigure(); plot(1e6*result.times,real(result.value),'b-',...
                1e6*result.times,real(exact.value),'r--');
xlabel('total echo time / microseconds'); ylabel('normalised coherence');
legend('64 sampled states','exhaustive product density'); kgrid();

end

