% Tests once-only thermalisation of additive SRSK relaxation. Syntax:
%
%                    result=test_srsk_thermal()
%
% Outputs:
%
%     result  - regression checks for rates, retention, and equilibrium
%
% A rapidly relaxing 14N source is coupled to 1H with positive, negative,
% or zero scalar coupling. An oriented quadrupole supplies a complex,
% noncommuting thermalisation Hamiltonian as a contrary control.
%
% talos@spindynamics.org

function result=test_srsk_thermal()

% State the regression target
result=new_test_result('kernel/srsk_thermal','SRSK thermalisation',...
                      'SRSK must add zero-destination rates before final thermalisation.');

% Specify a fast quadrupolar source and its proton partner
sys.magnet=1e-3; sys.isotopes={'1H','14N'};
inter.zeeman.scalar={0,0}; inter.coupling.scalar=cell(2);
inter.coupling.matrix=cell(2);
inter.coupling.matrix{2,2}=diag([-700 -300 1000]);
inter.relaxation={'t1_t2','SRSK'}; inter.srsk_sources=2;
inter.r1_rates={1,1e5}; inter.r2_rates={1,1e5};
inter.equilibrium='zero'; inter.temperature=1;
inter.rlx_keep='labframe'; inter.rlx_dfs='keep';
bas.formalism='sphten-liouv'; bas.approximation='none';
keep_modes={'labframe','diagonal','kite','secular'};
methods={'IME','dibari'}; angles=[0.3 0.7 0.2];

% Include the uncoupled limit and a sign reversal of the scalar coupling
for coupling=[100 -100 0]
    inter.coupling.scalar{1,2}=coupling;
    spin_system=test_spin_system(sys,inter,bas);
    unit=unit_state(spin_system);
    [H,Q]=hamiltonian(assume(spin_system,'labframe'),'left');

    % Test the isotropic call without an orientation argument
    R=relaxation(spin_system); rho_eq=equilibrium(spin_system,H);
    for n=1:numel(methods)
        spin_system.rlx.equilibrium=methods{n};
        if strcmp(methods{n},'IME')
            ref_thermal=thermalize(spin_system,R,[],[],rho_eq,'IME');
        else
            ref_thermal=thermalize(spin_system,R,H,1,[],'dibari');
        end
        observed=relaxation(spin_system);
        result=test_close(result,['isotropic ' methods{n}],observed,ref_thermal,1e-8,1e-12,...
                          'the isotropic call must also thermalise the total exactly once');
    end
    spin_system.rlx.equilibrium='zero';

    % Add the oriented quadrupolar interaction to the thermalisation target
    H=H+orientation(Q,angles); rho_eq=equilibrium(spin_system,H);
    assert(norm(imag(H),'fro')>1);

    % Independently augment the rates using the scalar-coupling formula
    delta_omega=diff(spin_system.inter.basefrqs);
    rate_long=(4/3)*(2*pi*coupling)^2*1e-5/(1+delta_omega^2*1e-10);
    rate_trans=(2/3)*(2*pi*coupling)^2*1e-5+rate_long/2;
    reference=spin_system; reference.rlx.theories={'t1_t2'};
    reference.rlx.r1_rates={1+rate_long,1e5};
    reference.rlx.r2_rates={1+rate_trans,1e5};

    % Preserve each supported spherical-tensor retention policy
    for k=1:numel(keep_modes)
        spin_system.rlx.keep=keep_modes{k};
        reference.rlx.keep=keep_modes{k};
        spin_system.rlx.equilibrium='zero';
        R=relaxation(spin_system,angles);
        ref_rates=relaxation(reference,angles);
        label=[keep_modes{k} ', J=' num2str(coupling)];
        result=test_close(result,['rates ' label],R,ref_rates,1e-8,1e-12,...
                          'zero-destination SRSK must preserve the physical source and added rates');
        assert(norm(R*H-H*R,'fro')>1);

        % Compare both thermalisation methods to the once-only reference
        for n=1:numel(methods)
            spin_system.rlx.equilibrium=methods{n};
            if strcmp(methods{n},'IME')
                ref_thermal=thermalize(spin_system,R,[],[],rho_eq,'IME');
            else
                ref_thermal=thermalize(spin_system,R,H,1,[],'dibari');
            end
            observed=relaxation(spin_system,angles);
            reference.rlx.equilibrium=methods{n};
            ref_rates=relaxation(reference,angles);
            reference.rlx.equilibrium='zero';
            result=test_close(result,['no SRSK ' label],observed,ref_rates,1e-8,1e-12,...
                              'explicit augmented rates without SRSK must give the same thermalisation');
            result=test_close(result,[methods{n} ' ' label],observed,ref_thermal,1e-8,1e-12,...
                              'the total generator must be thermalised exactly once');
            result=test_close(result,['equilibrium ' label],observed*rho_eq,0*rho_eq,1e-8,0,...
                              'the laboratory-frame thermal equilibrium must be stationary');
            result=test_close(result,['trace ' label],unit'*observed,0*unit',1e-8,0,...
                              'thermalisation must preserve the trace functional');
        end
    end
end

end


