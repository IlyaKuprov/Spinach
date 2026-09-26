% Tests SRSK formalism restrictions and supported source relaxation. Syntax:
%
%                         result=test_srsk_form()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% A rapidly relaxing 14N source broadens its scalar-coupled proton in
% spherical-tensor Liouville space. Zeeman Liouville SRSK must be refused
% explicitly rather than selecting a different high-spin relaxation model.
%
% talos@spindynamics.org

function result=test_srsk_form()

% State the formalism contract
result=new_test_result('kernel/srsk_form','SRSK formalism restriction',...
                       'SRSK requires the spherical-tensor Liouville formalism.');

% Build an applicable rapidly relaxing nitrogen source and proton partner
sys.magnet=1e-3; sys.isotopes={'1H','14N'};
inter.zeeman.scalar={0,0};
inter.coupling.scalar=cell(2); inter.coupling.scalar{1,2}=100;
inter.relaxation={'lindblad','SRSK'}; inter.srsk_sources=2;
inter.lind_r1_rates=[1,1e5]; inter.lind_r2_rates=[1,1e5];
inter.equilibrium='zero'; inter.temperature=298;
inter.rlx_keep='labframe'; inter.rlx_dfs='keep';
bas.formalism='zeeman-liouv'; bas.approximation='none';
spin_z=test_spin_system(sys,inter,bas);

% Require the explicit refusal across retention and equilibrium choices
retentions={'labframe','secular','diagonal'};
equilibria={'zero','IME','dibari'};
for n=1:numel(retentions)
    for k=1:numel(equilibria)
        spin_z.rlx.keep=retentions{n};
        spin_z.rlx.equilibrium=equilibria{k};
        message='';
        try
            relaxation(spin_z);
        catch err
            message=err.message;
        end
        result=test_true(result,['Zeeman ' retentions{n} ' ' equilibria{k}],...
                         contains(message,'SRSK')&&contains(message,'not implemented')&&...
                         contains(message,'zeeman-liouv'),...
                         'unsupported SRSK must be refused before constructing relaxation terms');
    end
end

% Keep the SRSK refusal specific even when extended T1/T2 is requested
spin_z.rlx.theories={'t1_t2','SRSK'};
message='';
try
    relaxation(spin_z);
catch err
    message=err.message;
end
result=test_true(result,'Zeeman T1/T2 plus SRSK',...
                 contains(message,'SRSK')&&contains(message,'not implemented'),...
                 'the public SRSK guard must precede the recursive model restriction');

% Preserve the supported Zeeman Lindblad control without SRSK
spin_z.rlx.theories={'lindblad'};
spin_z.rlx.keep='labframe'; spin_z.rlx.equilibrium='zero';
R=relaxation(spin_z);
rho=state(spin_z,'L+','1H');
result=test_close(result,'Zeeman transverse decay',R*rho,-rho,1e-10,1e-10,...
                  'removing SRSK leaves the specified proton transverse rate unchanged');

% Obtain the additive SRSK contribution in the supported formalism
bas.formalism='sphten-liouv';
spin_s=test_spin_system(sys,inter,bas);
R=relaxation(spin_s);
spin_s.rlx.theories={'lindblad'};
R=R-relaxation(spin_s);

% Check Abragam rates on longitudinal and complex non-Hermitian states
freq_diff=spin_s.inter.basefrqs(1)-spin_s.inter.basefrqs(2);
coupling=2*pi*100;
r1_add=(4/3)*coupling^2*(1e-5/(1+freq_diff^2*1e-10));
r2_add=(2/3)*coupling^2*(1e-5+1e-5/(1+freq_diff^2*1e-10));
rho_z=state(spin_s,'Lz','1H');
rho_p=(1+2i)*state(spin_s,'L+','1H');
result=test_close(result,'Spherical longitudinal SRSK',R*rho_z,-r1_add*rho_z,...
                  1e-10,1e-10,'the supported nitrogen source adds the Abragam longitudinal rate');
result=test_close(result,'Spherical transverse SRSK',R*rho_p,-r2_add*rho_p,...
                  1e-10,1e-10,'the supported nitrogen source adds the Abragam transverse rate');

end


