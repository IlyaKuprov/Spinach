% Tests diagonal relaxation retention formalism boundaries. Syntax:
%
%                    result=test_diagonal_guard()
%
% Outputs:
%
%     result  - regression result with explanatory messages
%
% Zeeman diagonal retention is not implemented; spherical-tensor diagonal
% and Zeeman full retention must preserve trace and the specified rates.
% A complex coherent state and noncommuting Hamiltonian exercise the full
% generator rather than relying only on population or eigenstate decay.
%
% talos@spindynamics.org

function result=test_diagonal_guard()

% State the retention boundary under test
result=new_test_result('kernel/diagonal_guard',...
                       'Diagonal relaxation retention guard',...
                       'unsupported Zeeman diagonal retention must fail without changing supported generators.');

% Specify spin-half and spin-one Lindblad baths
isotopes={'1H','14N'};
inter.zeeman.scalar={0}; inter.relaxation={'lindblad'};
inter.lind_r1_rates=4; inter.lind_r2_rates=7;
inter.equilibrium='zero'; inter.temperature=298;
inter.rlx_keep='labframe'; inter.rlx_dfs='keep';
sys.magnet=14.1;
bas.approximation='none';
for n=1:numel(isotopes)

    % Build the full Zeeman generator and its normalised identity
    sys.isotopes=isotopes(n); bas.formalism='zeeman-liouv';
    spin_system=test_spin_system(sys,inter,bas);
    R=relaxation(spin_system); unit=unit_state(spin_system);
    result=test_close(result,'Zeeman trace conservation',unit'*R,0*unit',1e-12,1e-12,...
                      'full retention preserves the left-null trace functional');
    result=test_close(result,'Zeeman identity stationarity',R*unit,0*unit,1e-12,1e-12,...
                      'the symmetric Lindblad bath leaves identity stationary');

    % Preserve the existing labframe default when retention is omitted
    spin_default=spin_system; spin_default.rlx=rmfield(spin_default.rlx,'keep');
    result=test_close(result,'omitted retention',relaxation(spin_default),R,0,0,...
                      'omitting retention must still select the full generator');

    % Propagate a complex pure state under noncommuting coherent dynamics
    bas_hilb.formalism='zeeman-hilb'; bas_hilb.approximation='none';
    spin_hilb=basis(spin_system,bas_hilb);
    H=0.3*operator(spin_hilb,'Lx',1)+0.2*operator(spin_hilb,'Ly',1);
    ket=zeros(size(H,1),1); ket(1)=1/sqrt(2); ket(2)=1i/sqrt(2);
    rho=ket*ket'; Q=hilb2liouv(H,'comm');
    result=test_true(result,'noncommuting fixture',norm(Q*R-R*Q,'fro')>1e-3,...
                     'coherent and dissipative dynamics must not commute in this control');
    rho_final=reshape(expm(full(-1i*Q+R))*rho(:),size(H));
    result=test_close(result,'coherent trace',trace(rho_final),1,1e-12,1e-12,...
                      'full retention preserves trace during complex coherent evolution');
    result=test_close(result,'coherent Hermiticity',rho_final,rho_final',1e-12,1e-12,...
                      'the supported full generator preserves Hermiticity');

    % Require the explicit unsupported-combination error
    spin_system.rlx.keep='diagonal'; message='';
    try
        relaxation(spin_system);
    catch err
        message=err.message;
    end
    result=test_true(result,'Zeeman diagonal refusal',...
                     contains(message,'not implemented')&&contains(message,'zeeman-liouv'),...
                     'both spin-half and spin-one diagonal Zeeman requests must be rejected explicitly');

    % Retain spherical-tensor self-relaxation and its physical rates
    bas.formalism='sphten-liouv'; spin_system=test_spin_system(sys,inter,bas);
    spin_system.rlx.keep='diagonal'; R=relaxation(spin_system);
    unit=unit_state(spin_system); rho_z=state(spin_system,'Lz',1);
    rho_p=state(spin_system,'L+',1);
    result=test_close(result,'spherical trace conservation',unit'*R,0*unit',1e-12,1e-12,...
                      'spherical-tensor diagonal retention preserves trace');
    result=test_close(result,'spherical identity stationarity',R*unit,0*unit,1e-12,1e-12,...
                      'spherical-tensor diagonal retention protects the unit state');
    result=test_close(result,'spherical longitudinal rate',R*rho_z,-4*rho_z,1e-12,1e-12,...
                      'the supported diagonal path preserves the specified R1 rate');
    result=test_close(result,'spherical transverse rate',R*rho_p,-7*rho_p,1e-12,1e-12,...
                      'the supported diagonal path preserves the specified R2 rate');
end

end


