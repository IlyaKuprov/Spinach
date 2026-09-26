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
% A shipped GISSMO XML fixture checks pure damping in both Liouville bases
% and the conversion from Lorentzian linewidth to transverse decay rate.
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

    % Preserve non-selective damping in both supported Liouville formalisms
    inter_damp=rmfield(inter,{'lind_r1_rates','lind_r2_rates'});
    inter_damp.relaxation={'damp'}; inter_damp.damp_rate=5;
    formalisms={'sphten-liouv','zeeman-liouv'};
    for k=1:numel(formalisms)
        bas.formalism=formalisms{k};
        spin_system=test_spin_system(sys,inter_damp,bas);
        R=relaxation(spin_system); unit=unit_state(spin_system);
        R_ref=-inter_damp.damp_rate*(unit_oper(spin_system)-unit*unit');
        R_ref=clean_up(spin_system,R_ref,spin_system.tols.rlx_zero);
        result=test_close(result,'full-retention damping',R,R_ref,1e-12,1e-12,...
                          'damp-only full retention preserves non-selective decay and the unit state');
    end
end

% Import an actual GISSMO subsystem with a one-hertz Lorentzian linewidth
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
[sys,inter]=gissmo2spinach(fullfile(root_dir,'examples','nmr_metabol','molecule_b.xml'),2);
result=test_close(result,'GISSMO linewidth rate',inter.damp_rate,pi,0,0,...
                  'the shipped one-hertz FWHM must give a pi-per-second damping rate');
result=test_true(result,'GISSMO full retention',strcmp(inter.rlx_keep,'labframe'),...
                 'pure damping must not request unsupported Zeeman diagonal retention');
bas.approximation='none'; formalisms={'sphten-liouv','zeeman-liouv'};
for n=1:numel(formalisms)

    % Verify the entire imported generator in both full Liouville bases
    bas.formalism=formalisms{n}; spin_system=test_spin_system(sys,inter,bas);
    R=relaxation(spin_system); unit=unit_state(spin_system);
    R_ref=-pi*(unit_oper(spin_system)-unit*unit');
    round_bound=spin_system.tols.rlx_zero*sqrt(nnz(R_ref))/2;
    R_ref=clean_up(spin_system,R_ref,spin_system.tols.rlx_zero);
    result=test_close(result,'GISSMO damping generator',norm(R-R_ref,'fro'),0,1e-12,1e-12,...
                      'the XML linewidth must damp every traceless state at the specified rate');
    result=test_close(result,'GISSMO trace conservation',unit'*R,0*unit',round_bound,1e-12,...
                      'the imported damping generator must conserve trace within entrywise rounding');
    result=test_close(result,'GISSMO identity stationarity',R*unit,0*unit,round_bound,1e-12,...
                      'the imported identity must remain stationary within entrywise rounding');

    % Recover the linewidth from the actual transverse decay eigenvalue
    rho_p=state(spin_system,'L+','1H');
    rate=-real(rho_p'*R*rho_p)/(rho_p'*rho_p);
    result=test_close(result,'GISSMO transverse decay',R*rho_p,-pi*rho_p,round_bound*norm(rho_p),1e-12,...
                      'the transverse signal decays as exp(-pi*FWHM*t)');
    result=test_close(result,'GISSMO linewidth in hertz',rate/pi,1,round_bound/pi,1e-12,...
                      'Lorentzian FWHM is the transverse decay rate divided by pi');

    % Preserve the formerly supported spherical diagonal-retention result
    if strcmp(bas.formalism,'sphten-liouv')
        spin_system.rlx.keep='diagonal'; R_diag=relaxation(spin_system);
        result=test_close(result,'GISSMO spherical preservation',norm(R-R_diag,'fro'),0,0,0,...
                          'damping is added after retention, so the spherical result must not change');
    end
end

end


