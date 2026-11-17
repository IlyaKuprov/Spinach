% Tests saturated CW detection against Bloch equations and a native
% propagator steady state, including concentrations, spatial embedding,
% the iterative branch, and rejected inputs. Syntax:
%
%                     result=test_slowpass_sat()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_slowpass_sat()

% Test record
result=new_test_result('kernel/slowpass_sat','Saturated CW steady states',...
                       'Checks Bloch physics, conserved concentrations, spatial embedding, and solver paths.');

% Isotropic electron with independently specified relaxation rates
sys.magnet=0.336; sys.isotopes={'E'};
inter.relaxation={'lindblad'}; inter.rlx_keep='labframe';
inter.lind_r1_rates=2e5; inter.lind_r2_rates=3e5;
inter.equilibrium='IME'; inter.temperature=298;
bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'esr');
H=hamiltonian(spin_system); R=relaxation(spin_system); K=kinetics(spin_system);
Lx=operator(spin_system,'Lx','E'); Lz=operator(spin_system,'Lz','E');
H=H+2*pi*1.2e5*Lz;
parameters.spins={'E'}; parameters.coil=coil_state(spin_system,'L+','E');
parameters.sweep=[-2e6 2e6]; parameters.npoints=17;
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));
z_eq=real(coil_state(spin_system,'Lz','E')'*rho_eq);
freq=2*pi*linspace(parameters.sweep(1),parameters.sweep(2),parameters.npoints)';

% Compare both quadratures across the saturation range
for power=[0 1e2 1e5 1e6]

    % Analytic Bloch response for the actual detuning and equilibrium
    parameters.irr_powers=power; omega=2*pi*power; delta=freq+2*pi*1.2e5;
    actual=slowpass_sat(spin_system,parameters,H,R,K);
    my=-omega*inter.lind_r2_rates*z_eq./...
       (inter.lind_r2_rates^2+delta.^2+omega^2*inter.lind_r2_rates/inter.lind_r1_rates);
    mx=-delta.*my/inter.lind_r2_rates;
    result=test_close(result,['Bloch nu1=' num2str(power)],actual,mx-1i*my,1e-12,1e-10,...
                      'both quadratures must obey the driven Bloch equations');

end

% Check an independent propagator fixed point at every offset
reference=zeros(size(freq));
for n=1:numel(freq)

    % Native steady state of the full thermalised generator
    L=H+1i*R+1i*K+omega*Lx+freq(n)*Lz;
    P=expm(full(-1i*L*5e-5));
    rho=steady(spin_system,P,[],'newton');
    reference(n)=parameters.coil'*rho;
    result=test_true(result,'stationary residual',norm(L*rho)<1e-7,...
                     'the independent fixed point must solve the generator equation');

end
result=test_close(result,'native steady state',actual,reference,1e-12,1e-10,...
                  'generator and propagator methods must agree');

% Retain the known trace contribution in a general observable
saved_coil=parameters.coil; parameters.coil=saved_coil+unit_state(spin_system);
with_unit=slowpass_sat(spin_system,parameters,H,R,K);
result=test_close(result,'identity component',with_unit,actual+1,1e-12,1e-10,...
                  'pinning must not remove the observable trace contribution');
parameters.coil=saved_coil;

% An initial state is deliberately irrelevant to the steady-state source
parameters.rho0=zeros(size(parameters.coil));
result=test_close(result,'rho0 ignored',slowpass_sat(spin_system,parameters,H,R,K),actual,1e-12,1e-10,...
                  'thermalised relaxation supplies the source rather than rho0');
parameters=rmfield(parameters,'rho0');

% Exercise GMRES above the preserved active-dimension threshold
spc_dim=25001; space_unit=sparse(1,1,1,spc_dim,1);
params_iter=parameters; params_iter.spc_dim=spc_dim; params_iter.npoints=3;
params_iter.coil=kron(space_unit,parameters.coil);
H_big=kron(speye(spc_dim),H);
R_big=kron(speye(spc_dim),R)-kron(spdiags([0;1e7*ones(spc_dim-1,1)],0,spc_dim,spc_dim),speye(4));
K_big=sparse(size(H_big,1),size(H_big,2));
iterated=slowpass_sat(spin_system,params_iter,H_big,R_big,K_big);
sample_idx=round(linspace(1,numel(actual),params_iter.npoints));
result=test_close(result,'iterative spatial solve',iterated,actual(sample_idx),1e-11,1e-9,...
                  'the iterative branch must preserve the independent rank-zero response');
clear H_big R_big K_big;

% Compare native liquid and gridfree isotropic limits with anisotropic diffusion
params_sle=parameters; params_sle.max_rank=2; params_sle.tau_c=diag([1e-8 2e-8 3e-8]);
params_sle.offset=0; params_sle.decouple={};
liquid_spec=liquid(spin_system,@slowpass_sat,parameters,'esr');
sle_spec=gridfree(spin_system,@slowpass_sat,params_sle,'esr');
result=test_close(result,'liquid and gridfree',sle_spec,liquid_spec,1e-12,1e-10,...
                  'free rotational diffusion must not alter an isotropic spin spectrum');
spin_zte=spin_system; spin_zte.sys.enable=[spin_zte.sys.enable {'zte'}];
result=test_close(result,'ZTE not invoked',gridfree(spin_zte,@slowpass_sat,params_sle,'esr'),sle_spec,1e-12,1e-10,...
                  'the saturated solver does not call trajectory reduction');

% Exercise finite-input checks without accepting a downstream numerical failure
bad_values={NaN,Inf,-1,[1 2]};
for n=1:numel(bad_values)

    % Reject a malformed irradiation amplitude in the grumbler
    bad=parameters; bad.irr_powers=bad_values{n}; rejected=false;
    try
        slowpass_sat(spin_system,bad,H,R,K);
    catch err
        rejected=contains(err.message,'parameters.irr_powers');
    end
    result=test_true(result,'invalid nutation frequency',rejected,'invalid driving must be rejected');

end
bad=parameters; bad.sweep=[0 NaN]; rejected=false;
try
    slowpass_sat(spin_system,bad,H,R,K);
catch err
    rejected=contains(err.message,'parameters.sweep');
end
result=test_true(result,'nonfinite sweep',rejected,'nonfinite offsets must be rejected');

% Reject zero-equilibrium relaxation and population-changing kinetics
R_zero=R; R_zero(:,1)=0; rejected=false;
try
    slowpass_sat(spin_system,parameters,H,R_zero,K);
catch err
    rejected=contains(err.message,'not thermalised');
end
result=test_true(result,'unthermalised relaxation',rejected,'a thermal source is required');
K_loss=K-speye(size(K)); rejected=false;
try
    slowpass_sat(spin_system,parameters,H,R,K_loss);
catch err
    rejected=contains(err.message,'not conserved');
end
result=test_true(result,'population loss',rejected,'individual substance populations must be conserved');

% Verify concentration weighting for two independent substances
single_spec=slowpass_sat(spin_system,parameters,H,R,K);
sys.isotopes={'E','E'}; inter.chem.parts={1,2}; inter.chem.concs=[0.2 1.7];
inter.lind_r1_rates=[2e5 2e5]; inter.lind_r2_rates=[3e5 3e5];
bas.approximation={'none','none'};
mixture=test_spin_system(sys,inter,bas); mixture=assume(mixture,'esr');
H_mix=hamiltonian(mixture)+2*pi*1.2e5*operator(mixture,'Lz','E');
R_mix=relaxation(mixture); K_mix=kinetics(mixture);
params_mix=parameters; params_mix.coil=coil_state(mixture,'L+','E');
mixed=slowpass_sat(mixture,params_mix,H_mix,R_mix,K_mix);
result=test_close(result,'unequal concentrations',mixed,1.9*single_spec,1e-12,1e-10,...
                  'independent substance signals must be weighted once');
mixture.chem.concs=[0 1.7];
result=test_close(result,'zero concentration',slowpass_sat(mixture,params_mix,H_mix,R_mix,K_mix),...
                  1.7*single_spec,1e-12,1e-10,'zero population must produce no signal');

end


