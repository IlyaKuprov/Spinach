% Focused checks of the ZTE leakage-based time estimate and pruning.
%
% Syntax: result=test_zte_warranty()
%
% Output: regression result; assertions fail on an incorrect estimate.
%
% ilya.kuprov@weizmann.ac.il

function result=test_zte_warranty()

% Set independent selection and estimate tolerances
spin_system.sys.output='hush'; spin_system.sys.enable={}; spin_system.sys.disable={};
spin_system.bas.formalism='sphten-liouv';
[spin_system,~]=tolerances(spin_system,struct());
spin_system.tols.zte_nsteps=2; spin_system.tols.zte_maxden=1;
spin_system.tols.zte_tol=1e-3;

% Check the linear leakage estimate and its small-time unitary reference
L=sparse([0 1e-5;1e-5 0]); rho=[1;0];
P=zte(spin_system,L,rho,1);
duration=zte_warr(spin_system,L,rho,P);
assert(abs(duration-0.1)<eps);
assert(norm(expm(-1i*full(L)*duration)*rho-P*expm(-1i*full(P'*L*P)*duration)*(P'*rho))<=1e-6);
assert(abs(zte_warr(spin_system,L,[1;5e-7],P)-0.05)<eps);

% Check unchanged, invariant, empty, scalar, and initial-discard cases
assert(isinf(zte_warr(spin_system,L,rho,1)));
assert(isinf(zte_warr(spin_system,L,rho,speye(2))));
assert(isinf(zte_warr(spin_system,speye(2),rho,P)));
assert(isinf(zte_warr(spin_system,L,[0;0],sparse(2,0))));
assert(isinf(zte_warr(spin_system,2,1,1)));
assert(zte_warr(spin_system,L,[1;1e-6],P)==0);
assert(zte_warr(spin_system,L,[1;1e-5],P)==0);
assert(zte_warr(spin_system,L,rho,sparse(2,0))==0);

% Check zero-state, density, disabled, and small-state skip paths
assert(isequal(zte(spin_system,L,[0;0]),1));
spin_system.tols.zte_maxden=0.5;
assert(isequal(zte(spin_system,L,[1;1]),1));
spin_system.sys.disable={'zte'};
assert(isequal(zte(spin_system,L,rho),1));
spin_system.sys.disable={}; spin_system.tols.zte_maxden=1;
assert(isequal(zte(spin_system,L,[1e-8;1e-8]),1));
assert(isequal(zte(spin_system,L,rho,2),speye(2)));

% Verify that the estimate tolerance never changes selection
P=zte(spin_system,L,rho);
spin_system.tols.zte_warr=1e-12;
assert(isequal(zte(spin_system,L,rho),P));
spin_system.tols.zte_warr=1;
assert(isequal(zte(spin_system,L,rho),P));
spin_system.tols.zte_warr=1e-6;

% Check the explicitly labelled diagnostic
spin_system.sys.output=1;
text=evalc('zte(spin_system,L,rho);');
assert(contains(text,'ZTE warranty estimate: absolute 2-norm tolerance'));
assert(contains(text,'estimated time')&&contains(text,'not a guarantee'));

% Reject invalid warranty tolerances and accept the standard override route
spin_system.sys.output='hush';
invalid={0,-1,NaN,Inf,1i,[1 2],'bad',true};
for n=1:numel(invalid)
    sys.tols.zte_warr=invalid{n}; failed=false;
    try
        tolerances(spin_system,sys);
    catch exception
        failed=contains(exception.message,'zte_warr');
    end
    assert(failed);
end
sys.tols.zte_warr=2e-4;
[updated,parsed]=tolerances(spin_system,sys);
assert(updated.tols.zte_warr==2e-4&&~isfield(parsed,'tols'));

% Return the regression result after all assertions pass
result=new_test_result('kernel/zte_warranty','ZTE supplied-vector warranty',...
                       'Leakage time estimates, skipped reductions, and tolerance independence.');

end


