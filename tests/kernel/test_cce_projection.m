% Test: kernel/cce_projection | CCE observed-spin projection
% Tests observed expectations and virtual-level corrections. Syntax:
%
%                    result=test_cce_projection()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_projection()

% State the independent projection targets
result=new_test_result('kernel/cce_projection','CCE observed projection',...
                       'Check explicit states, projected fields and resolved virtual levels.');

% Create an observed spin-half with an exactly known longitudinal gap
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
system=create(sys,struct());
gap=1000;
coupling=4;
system.inter.zeeman.matrix{1}=gap*eye(3);
system=cce_test_write(system,[1 2],diag([coupling 0 0]));
parameters.observed=1;
parameters.cce.bath_spins=2;
parameters.times=0;
parameters.orientation=[0 0 0];
parameters.cce.interaction_cutoff=Inf;
parameters.cce.numerics=struct('atol',1e-10,'rtol',1e-12,'max_memory_bytes',1e8);
model=cce_compile(system,parameters,'labframe');
projection=struct('method','cce','transition',struct('alpha',[1;0],'beta',[0;1]));
manifolds=cce_project(model,projection);
result=test_close(result,'alpha expectation',...
                  manifolds.alpha.expectations,[0 0 1/2],0,0,'Spin-half up state has positive one-half Sz.');
result=test_close(result,'beta expectation',...
                  manifolds.beta.expectations,[0 0 -1/2],0,0,'Spin-half down state has negative one-half Sz.');
result=test_true(result,'second order off by default',...
                 ~manifolds.second_order_enabled&&norm(manifolds.alpha.susceptibility,'fro')==0,'No implicit perturbative correction.');

% Compare manifold Hamiltonians against explicit full-Hilbert projections
[H,meta]=cce_hamiltonian(model,2,struct());
for name={'alpha','beta'}
    manifold=manifolds.(name{1});
    embedding=kron(manifold.vector,eye(2));
    projected=embedding'*H*embedding;
    reference=manifold.energy*eye(2)+coupling*manifold.expectations(1)*[0 1;1 0]/2;
    result=test_close(result,[name{1} ' projected projection'],...
                      projected,reference,meta.roundoff_bound,0,'Independent dense product-space projection.');
end

% Enable an observed-gap expansion and compare against the exact coupled spectrum
projection.second_order=struct('enabled',true,'gap_tolerance',1,'overlap_tolerance',1e-12);
manifolds=cce_project(model,projection);
result=test_close(result,'virtual coefficient',...
                  manifolds.alpha.susceptibility(1,1),1/(4*gap),eps,0,'Analytic spin-half transverse susceptibility.');
correction=coupling^2*manifolds.alpha.susceptibility(1,1)/4;
exact=sort(eig(full(H)));
approximate=gap/2+correction;
error_value=abs(approximate-exact(end));
bound=coupling^4/(256*gap^3);
result=test_true(result,'second-order spectrum',...
                 error_value<=bound+eps(gap),'Exact square-root spectrum bounds omitted fourth order.');
fprintf('CCE projection: exact spectral error=%g, fourth-order bound=%g\n',error_value,bound);

% Reject a near-degenerate virtual denominator before any division
projection.second_order.gap_tolerance=2*gap;
rejected=false;
try
    cce_project(model,projection);
catch
    rejected=true;
end
result=test_true(result,'near-degenerate denominator rejected',...
                 rejected,'An explicit gap guard prevents unresolved perturbation theory.');

% Reject non-eigenstate projectors and gCCE double counting
projection.second_order.gap_tolerance=1;
projection.transition.alpha=[1;1]/sqrt(2);
projection.transition.beta=[1;-1]/sqrt(2);
rejected=false;
try
    cce_project(model,projection);
catch
    rejected=true;
end
result=test_true(result,'overlap guard',...
                 rejected,'Explicit superposition is not silently reselected.');
projection.transition.alpha=[1;0];
projection.transition.beta=[0;1];
projection.method='gcce';
rejected=false;
try
    cce_project(model,projection);
catch
    rejected=true;
end
result=test_true(result,'gCCE correction rejected',...
                 rejected,'Full observed dynamics already contains virtual transitions.');

% Verify multi-centre product ordering and a complex orthogonal selected two-level transition basis
projection=rmfield(projection,'second_order');
parameters.observed=[2 1];
parameters.cce.bath_spins=[];
projection.transition.alpha=kron([1;0],[1;1i]/sqrt(2));
projection.transition.beta=kron([0;1],[1;-1i]/sqrt(2));
manifolds=cce_project(cce_compile(system,parameters,'labframe'),projection);
result=test_close(result,'multiple observed-spin system',...
                  manifolds.alpha.expectations,[0 0 1/2;0 1/2 0],8*eps,0,'Complex state expectations preserve declared observed order.');

% Compare mediated pair dynamics against independently projected dense matrices
sys.isotopes={'1H','1H','1H'};
system=create(sys,struct());
system.inter.zeeman.matrix{1}=1000*eye(3);
system=cce_test_write(system,[1 2],diag([100 0 0]));
system=cce_test_write(system,[1 3],diag([80 0 0]));
p.observed=1;
p.times=[0 0.1 0.3 0.5];
p.orientation=[0 0 0];
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.cce=struct('bath_spins',[2 3],'approximation','projected',...
             'expansion','standard','order',2,'interaction_cutoff',Inf,...
             'meanfield','none');
p.cce.graph=struct('mode','coupling','threshold',0);
p.cce.bath_state=struct('mode','mixed_product');
p.cce.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',1e8);
p.cce.second_order=struct('enabled',true,'gap_tolerance',1,'overlap_tolerance',1e-10);
[actual,~]=clexpn(system,@free_evolution,p,'labframe');
x=[0 1;1 0]/2;
field=100*kron(x,eye(2))+80*kron(eye(2),x);
ha=500*eye(4)+field^2/4000;
hb=-ha;
reference=arrayfun(@(t) trace(expm(-1i*ha*t)*expm(1i*hb*t))/4,p.times).';
result=test_close(result,'second-order mediated graph',actual,reference,...
                  1e-10,0,'The coupling graph must retain susceptibility-induced bath pairs.');
fprintf('R7 projection sparse dense discrepancy %.17g\n',max(abs(actual-reference)));

% Preserve the explicit opt-out without altering the projected Hamiltonian
p.cce.graph.approximate=true;
[approximate,~]=clexpn(system,@free_evolution,p,'labframe');
result=test_true(result,'second-order graph approximation explicit',...
                 max(abs(approximate-reference))>0.4,...
                 'Only an explicitly approximate graph may omit the induced pair.');

% Compare the explicit midpoint Hahn approximation against two conditional histories
system.inter.zeeman.matrix{2}=630*eye(3);
system=cce_test_write(system,[1 2],[100 0 240; 0 0 0; 0 0 0]);
p.transition=struct('spins',1,'levels',[1 2]);
p.times=[0 0.001 0.003];
p.cce.bath_spins=2;
p.cce.order=1;
p.cce.second_order.enabled=false;
z=diag([1 -1])/2;
ha=630*z-120*x;
hb=630*z+120*x;
reference=ones(numel(p.times),1);
for k=2:numel(p.times)

    % Evaluate the two histories without invoking the native retained generator
    ua=expm(-0.5i*p.times(k)*ha);
    ub=expm(-0.5i*p.times(k)*hb);
    reference(k)=trace(ua*ub*(ub*ua)')/2;

end

% The same transverse coupling must still be rejected by invariant projection
rejected=false;
try
    clexpn(system,@hahn_echo,p,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce:population_transfer');
end
result=test_true(result,'invariant Hahn rejects register mixing',rejected,...
                 'The opt-in conditional approximation must not change the default.');

% Exercise both retained-Hamiltonian backends through the registered experiment
p.cce.projection='conditional';
for backend={'direct','stock'}

    % Check the complex response before the plotted real part or magnitude
    p.cce.backend=backend{1};
    [actual,details]=clexpn(system,@hahn_echo,p,'labframe');
    bound=norm(details.raw{1}(2).abs_error_bound)+numel(reference)*eps;
    result=test_close(result,['conditional Hahn ' backend{1}],actual,reference,bound,0,...
                      'The isolated register phase cancels but both bath terms remain.');

end

% An off-centre pulse would not cancel the isolated register phase
p.refocusing=struct('spins',1,'axis','x','angle',pi,'fraction',0.25);
rejected=false;
try
    clexpn(system,@hahn_echo,p,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce_options:projection');
end
result=test_true(result,'conditional Hahn rejects off-centre refocusing',rejected,...
                 'The approximation is restricted to the ideal midpoint pulse.');

end


