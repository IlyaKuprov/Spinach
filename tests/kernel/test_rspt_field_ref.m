% Tests Zeeman reference fields in perturbative eigensystems. Syntax:
%
%                    result=test_rspt_field_ref()
%
% Outputs:
%
%     result - field scaling, slopes, populations, and validation checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_rspt_field_ref()

% State the field-reference contract being tested
result=new_test_result('kernel/rspt_field_ref',...
                       'RSPT Zeeman reference fields',...
                       'Reference-field scaling must preserve energies, slopes, and populations.');

% Build a quiet two-level descriptor with a finite thermal scale
spin_system.bas.formalism='zeeman-hilb';
spin_system.bas.nsubst=1; spin_system.inter.magnet=1;
spin_system.sys.output='hush'; spin_system.sys.enable={};
spin_system.sys.disable={}; spin_system.chem.parts={1};
spin_system.chem.concs=1; spin_system.rlx.temperature=1;
spin_system.tols.hbar=1; spin_system.tols.kbol=1;
spin_system.tols.small_matrix=10; spin_system.tols.prop_chop=1e-12;
spin_system.tols.dense_matrix=0.5; spin_system.tols.liouv_zero=1e-12;
Hz=diag([-2 2]); Hc=[0.1 0.003;0.003 -0.1]; Hmw=[0 1;1 0];
parameters.rspt_order=Inf;

% Compare each reference field against the same direct Hamiltonian
for reference=[0.35 1 3.4]
    spin_system.inter.magnet=reference;
    field=0.7;
    [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,reference*Hz,Hc,Hmw,field);
    [V_ref,E_ref]=eig(field*Hz+Hc,'vector');
    [E_ref,idx]=sort(E_ref); V_ref=V_ref(:,idx);
    dE_ref=real(diag(V_ref'*Hz*V_ref));
    LP_ref=exp(-(E_ref-min(E_ref))); LP_ref=LP_ref/sum(LP_ref);
    result=test_close(result,'reference-field eigenenergies',E,E_ref,1e-12,1e-12,...
                      'Every Zeeman reference must produce H(B)=B*Hz+Hc.');
    result=test_close(result,'reference-field eigenvectors',abs(V),abs(V_ref),1e-12,1e-12,...
                      'Eigenvectors must agree up to their arbitrary signs.');
    result=test_close(result,'reference-field derivatives',dE,dE_ref,1e-12,1e-12,...
                      'Hellmann-Feynman slopes use the per-Tesla Zeeman operator.');
    result=test_close(result,'reference-field transitions',T,abs(V_ref'*Hmw*V_ref).^2,1e-12,1e-12,...
                      'Transition moments must not depend on the reference magnet.');
    result=test_close(result,'reference-field populations',LP,LP_ref,1e-12,1e-12,...
                      'Thermal populations must use the target-field Hamiltonian.');
end

% Check perturbative one-Tesla equivalence and convergence at other references
parameters.rho0=diag([0.8 0.2]);
for order=1:4
    parameters.rspt_order=order;
    spin_system.inter.magnet=1;
    [E_one,V_one,dE_one,T_one,LP_one]=...
        rspt_eig(spin_system,parameters,Hz,Hc,Hmw,field);
    [E_old,V_old]=rspert(diag(field*Hz+Hc),Hc-diag(diag(Hc)),order);
    [E_old,idx]=sort(E_old); V_old=V_old(:,idx);
    result=test_close(result,'one-Tesla perturbative contract',...
                      [E_one;abs(V_one(:))],[E_old;abs(V_old(:))],1e-12,1e-12,...
                      'At one Tesla, the original perturbative Hamiltonian is unchanged.');
    spin_system.inter.magnet=3.4;
    [E,~,dE,T,LP]=rspt_eig(spin_system,parameters,3.4*Hz,Hc,Hmw,field);
    result=test_close(result,'perturbative reference invariance',...
                      [E;dE;T(:);LP],[E_one;dE_one;T_one(:);LP_one],1e-12,1e-12,...
                      'All perturbative outputs must be independent of the reference field.');
    if order==4
        result=test_close(result,'fourth-order direct energies',E,E_ref,1e-11,1e-11,...
                          'Weakly mixed fourth-order energies must approach diagonalisation.');
    end
end

% Compare exact slopes with target-field finite differences
parameters.rspt_order=Inf;
step=1e-5;
[~,~,dE]=rspt_eig(spin_system,parameters,3.4*Hz,Hc,Hmw,field);
E_plus=rspt_eig(spin_system,parameters,3.4*Hz,Hc,Hmw,field+step);
E_minus=rspt_eig(spin_system,parameters,3.4*Hz,Hc,Hmw,field-step);
result=test_close(result,'target-field finite-difference slopes',...
                  dE,(E_plus-E_minus)/(2*step),1e-8,1e-8,...
                  'Reference scaling must also apply to derivatives away from the reference field.');

% Require undefined or invalid reference fields to fail explicitly
invalid={0,NaN,Inf,1i,[1 2],'one'};
for reference=invalid
    spin_system.inter.magnet=reference{1}; rejected=false;
    try
        rspt_eig(spin_system,parameters,Hz,Hc,Hmw,field);
    catch exception
        rejected=contains(exception.message,'finite non-zero real scalar');
    end
    result=test_true(result,'invalid reference rejection',rejected,...
                     'Invalid reference magnets must not silently scale the Zeeman term.');
end
spin_system.inter=rmfield(spin_system.inter,'magnet'); rejected=false;
try
    rspt_eig(spin_system,parameters,Hz,Hc,Hmw,field);
catch exception
    rejected=contains(exception.message,'finite non-zero real scalar');
end
result=test_true(result,'missing reference rejection',rejected,...
                 'The reference field is required rather than guessed.');

end


