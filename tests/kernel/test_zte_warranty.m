% ZTE warranty checks against dense exponential reference trajectories.
% Tests unitary, contractive, non-normal amplifying, asymmetric-exchange,
% invariant, and initial-discard cases, as well as skips and tolerances.
%
% Syntax: result=test_zte_warranty()
%
% Output: regression result; assertions fail on an incorrect warranty.
%
% ilya.kuprov@weizmann.ac.il

function result=test_zte_warranty()

% Set independent selection and warranty tolerances
spin_system.sys.output='hush'; spin_system.sys.enable={}; spin_system.sys.disable={};
spin_system.bas.formalism='sphten-liouv';
[spin_system,~]=tolerances(spin_system,struct());
spin_system.tols.zte_nsteps=2; spin_system.tols.zte_maxden=1;
spin_system.tols.zte_tol=1e-3; spin_system.sys.output=1;

% Cover unitary, contractive, amplifying, asymmetric, and invariant generators
matrices={[-1i 1e-5;-1e-5 -2i],[-1 0;1e-5 -2],...
          [2 1;1e-5 3],[-0.5 2;0.5 -2],diag([1 2]),zeros(2)};
states={[1;0],[1;0],[1;0],[1;0],[1;0],[1;1e-8]};
labels={'unitary','contractive','amplifying','asymmetric','invariant','initial_only'};
for n=1:numel(matrices)
    rho=states{n};
    text=evalc('P=zte(spin_system,1i*sparse(matrices{n}),rho,1);');
    token=regexp(text,'t <= ([^ ]+) seconds','tokens','once');
    duration=str2double(token{1});
    retained=find(any(P,2)); discarded=setdiff(1:2,retained);
    initial=norm(rho(discarded));
    leakage=norm(full(matrices{n}(discarded,retained)),2);
    growth=max(0,max(eig((matrices{n}+matrices{n}' )/2)));
    assert(initial<=spin_system.tols.zte_warr);
    if isinf(duration), times=linspace(0,1,21); else, times=linspace(0,duration,21); end
    errors=zeros(size(times));
    for k=1:numel(times)
        errors(k)=norm(expm(matrices{n}*times(k))*rho-P*expm(P'*matrices{n}*P*times(k))*(P'*rho));
    end
    assert(all(errors<=spin_system.tols.zte_warr+1e-13));
    fprintf('CASE %s warranty %.17g max_error %.17g numerical_abscissa %.9g leakage %.9g\n',...
            labels{n},duration,max(errors),growth,leakage);
end

% Check no interval, zero state, density skip, disabled ZTE, and no dropping
text=evalc('P=zte(spin_system,sparse(2,2),[1;1e-5],1);');
assert(contains(text,'no interval (not even 0 seconds)'));
text=evalc('P=zte(spin_system,sparse(2,2),[0;0]);');
assert(isequal(P,1)&&contains(text,'Inf seconds'));
spin_system.tols.zte_maxden=0.5;
text=evalc('P=zte(spin_system,speye(2),[1;1]);');
assert(isequal(P,1)&&contains(text,'Inf seconds'));
spin_system.sys.disable={'zte'};
text=evalc('P=zte(spin_system,speye(2),[1;0]);');
assert(isequal(P,1)&&contains(text,'Inf seconds'));
spin_system.sys.disable={}; spin_system.tols.zte_maxden=1;
text=evalc('P=zte(spin_system,speye(2),[1;1]);');
assert(isequal(P,speye(2))&&contains(text,'Inf seconds'));

% Exercise log-domain evaluation when a direct exponential would overflow
spin_system.tols.zte_nsteps=1;
text=evalc('P=zte(spin_system,1i*sparse([1e300 0;1e-300 1e300]),[1;0],1);');
token=regexp(text,'t <= ([^ ]+) seconds','tokens','once');
duration=str2double(token{1});
assert(isfinite(duration)&&(duration>0));
assert(1e300*duration+log(1e-300)+log(duration)<=log(1e-6)+1e-12);

% Check the zero-time boundary and the small-state skip
text=evalc('P=zte(spin_system,1i*speye(2),[1;1e-6],1);');
assert(contains(text,'t <= 0 seconds'));
spin_system.tols.zte_tol=1e-5;
text=evalc('P=zte(spin_system,1i*speye(2),[1e-8;1e-8]);');
assert(isequal(P,1)&&contains(text,'Inf seconds'));
spin_system.tols.zte_nsteps=2;

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
                       'Absolute state error stays within the reported time-domain budget.');

end


