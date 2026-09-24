% Tests reduction of true state columns with permutation symmetry. Syntax:
%
%                      result=test_stack_reduce()
%
% Outputs:
%
%     result - regression checks against independent propagation and
%              single-column screening references
%
% The tests cover antisymmetric, symmetric, and complex columns in three
% formalisms, dissipative propagation, ZTE ranking, and path tracing.
% Wide sparse and tolerance-relative mixed-scale stacks exercise a real
% parallel pool; non-unitary generators test growth and dynamical ranking.
%
% talos@spindynamics.org

function result=test_stack_reduce()

% Initialise the test result
result=new_test_result('kernel/stack_reduce','Stacked-state reduction',...
                       'screening must preserve every actual state column.');

% Build a physical pair of magnetically equivalent protons
sys.magnet=9.4; sys.isotopes={'1H','1H'};
sys.parallel={'processes',1}; sys.output='hush';
sys.disable={'hygiene'}; sys.parprops={};
inter.zeeman.scalar={0,0}; inter.coupling.scalar={0,10;10,0};
bas.approximation='none'; bas.sym_group={'S2'};
bas.sym_spins={[1 2]}; bas.sym_a1g_only=false;
forms={'sphten-liouv','zeeman-liouv','zeeman-wavef'};
for k=1:numel(forms)

    % Obtain the physical generator and opposite permutation sectors
    bas.formalism=forms{k};
    spin_system=basis(create(sys,inter),bas);
    spin_system=assume(spin_system,'nmr'); H=hamiltonian(spin_system);
    if strcmp(forms{k},'zeeman-wavef')
        rho=[0;1;-1;0]/sqrt(2); sym=[0;1;1;0]/sqrt(2);
    else
        rho=state(spin_system,'Lz',1)-state(spin_system,'Lz',2);
        sym=state(spin_system,'Lz',1)+state(spin_system,'Lz',2);
        rho=rho/norm(rho); sym=sym/norm(sym);
        H=H-0.3i*speye(size(H));
    end

    % Compare single, repeated, phase-shifted, and mixed-sector columns
    stacks={rho,[rho rho],[rho 1i*rho],[rho sym],...
            [rho+1i*sym rho-1i*sym zeros(size(rho))]};
    for n=1:numel(stacks)
        inputs=stacks{n}; exact=expm(full(-0.003i*H))*inputs;
        got=evolution(spin_system,H,[],inputs,0.003,1,'final');
        result=test_close(result,[forms{k} ' stack ' num2str(n)],got,exact,1e-7,1e-7,...
                          'every column must follow its own exact exponential trajectory');
    end

    % Require useful reduction rather than a blanket identity projector
    projectors=reduce(spin_system,H,[rho 1i*rho]);
    P=[projectors{:}];
    result=test_true(result,[forms{k} ' reduced dimension'],...
                     (size(P,1)==size(H,1))&&(size(P,2)<size(H,1)),...
                     'unoccupied symmetry sectors and coordinates must still be removed');
    result=test_close(result,[forms{k} ' retained columns'],P*(P'*[rho 1i*rho]),...
                      [rho 1i*rho],1e-8,1e-8,'projection must preserve both true columns');
end

% Build Liouville-space controls without permutation factorisation
bas=struct('formalism','sphten-liouv','approximation','none');
spin_system=basis(create(sys,inter),bas);
spin_system=assume(spin_system,'nmr'); H=hamiltonian(spin_system);
rho=state(spin_system,'L+',1)-state(spin_system,'L+',2);
sym=state(spin_system,'Lz',1)+state(spin_system,'Lz',2);
inputs=[rho 1i*rho sym];
spin_system.tols.zte_maxden=1;
P=zte(spin_system,H,inputs); retained=diag(P*P');
expected=false(size(H,1),1);
for n=1:size(inputs,2)

    % Form the union of separately screened actual trajectories
    Q=zte(spin_system,H,inputs(:,n));
    if isscalar(Q), Q=speye(size(H)); end
    expected=expected|logical(diag(Q*Q'));
end
result=test_close(result,'ZTE trajectory union',retained,double(expected),0,0,...
                  'stack screening must retain the union of individual trajectories');

% Check coordinate ranking and dimensions with a zero generator
inputs=sparse(size(H,1),size(H,1)+1);
inputs(2,1)=3; inputs(4,end)=2i; inputs(7,2)=-1;
H=sparse(size(H,1),size(H,1));
P=zte(spin_system,H,inputs);
Q=speye(size(H)); Q=Q(:,[2 4 7]);
result=test_close(result,'zero-generator support',P*P',Q*Q',0,0,...
                  'a wide stack must produce one coordinate mask, not a mask per column');
P=zte(spin_system,H,inputs,2);
Q=speye(size(H)); Q=Q(:,[2 4]);
result=test_close(result,'explicit state ranking',P*P',Q*Q',0,0,...
                  'nstates ranks maximum amplitude over columns and time');
P=zte(spin_system,H,inputs,1);
result=test_true(result,'one retained coordinate',size(P,2)==1,...
                 'the lower nstates boundary must retain exactly one coordinate');
P=zte(spin_system,H,inputs,size(H,1));
result=test_close(result,'full retained dimension',P*P',speye(size(H)),0,0,...
                  'the upper nstates boundary is the row dimension');

% Preserve the existing density and small-amplitude shortcuts
spin_system.tols.zte_maxden=0.5;
P=zte(spin_system,H,ones(size(inputs)),2);
result=test_true(result,'dense shortcut',isequal(P,1),...
                 'density screening retains its existing shortcut');
rho_scale=spin_system.tols.zte_tol/(2*norm(inputs,1));
P=zte(spin_system,H,inputs*rho_scale,2);
result=test_true(result,'tiny shortcut',isequal(P,1),...
                 'small-amplitude screening retains its existing shortcut');

% Reject a state count that confuses columns with basis coordinates
caught=false;
try
    zte(spin_system,H,inputs,size(H,1)+1);
catch exception
    caught=contains(exception.message,'state space dimension');
end
result=test_true(result,'nstates dimension guard',caught,...
                 'the number of columns must not enlarge the state space dimension');

% Preserve the explicit disable override
spin_system.sys.disable={'zte'};
P=zte(spin_system,H,inputs,2);
result=test_true(result,'ZTE disable',isequal(P,1),...
                 'an explicit disable must still leave the space unchanged');

% Force disconnected-subspace screening of complex columns
spin_system.sys.disable={'merge'};
spin_system.tols.merge_dim=1;
projectors=path_trace(spin_system,H,inputs);
P=[projectors{:}]; Q=speye(size(H)); Q=Q(:,[2 4 7]);
result=test_close(result,'path-tracing column union',P*P',Q*Q',0,0,...
                  'matrix 1-norm screening retains components occupied in any column');

% Exercise mixed scales under the real pool conditions of stack propagation
result=test_true(result,'active propagation pool',poolsize>0,...
                 'the wide-stack regression must run with an actual parallel pool');
weak=sqrt(spin_system.tols.zte_tol*eps('double'));
result=test_true(result,'tolerance-relative weak column',...
                 (weak>spin_system.tols.zte_tol)&&(weak<eps('double')),...
                 'the weak column lies above ZTE tolerance but below machine epsilon');
H=sparse([2 3],[3 2],1,300,300);
inputs=sparse([1 2 5],[1 2 3],[1 1i*weak spin_system.tols.zte_tol/2],300,40);
P=zte(spin_system,H,inputs);
expected=zeros(300,1); expected([1 2 3])=1;
result=test_close(result,'mixed-scale support',full(any(P,2)),expected,0,0,...
                  'weak populated columns must reach their coupled coordinates');
Q=zte(spin_system,H,full(inputs));
result=test_close(result,'dense mixed-scale support',full(any(Q,2)),expected,0,0,...
                  'screening must not depend on sparse versus full input storage');
control=spin_system; control.tols.zte_maxden=1;
one=zte(control,H,inputs(:,2));
result=test_true(result,'independent weak trajectory',...
                 isequal(find(any(one,2)),[2;3]),...
                 'the stack retains the same weak trajectory as independent screening');

% Retain nonzero columns that grow across the tolerance during propagation
H=sparse(2,2,1i,300,300);
inputs=sparse([1 2],[1 2],[1 spin_system.tols.zte_tol/2],300,40);
P=zte(spin_system,H,inputs);
result=test_true(result,'subthreshold initial column growth',...
                 isequal(find(any(P,2)),[1;2]),...
                 'only exactly zero columns may be omitted before propagation');

% Rank dynamical row maxima rather than norms or input-column amplitudes
H=sparse([2 3],[1 2],1,300,300);
inputs=sparse([1 4],[1 40],[1 0.8i],300,40);
P=zte(spin_system,H,inputs,2);
result=test_true(result,'dynamical state ranking',...
                 isequal(find(any(P,2)),[2;3]),...
                 'two nilpotent-generator steps give row maxima [1 2 2 0.8]');

% Screen a genuinely wide sparse stack without losing useful reduction
H=sparse([1 2],[2 1],1,4096,4096);
inputs=sparse(ones(1,512),1:512,ones(1,512),4096,513);
inputs(:,2:2:512)=1i*inputs(:,2:2:512);
P=zte(spin_system,H,inputs);
result=test_true(result,'wide sparse screening',...
                 issparse(inputs)&&isequal(size(P),[4096 2])&&...
                 isequal(find(any(P,2)),[1;2]),...
                 '512 phased columns and one zero column retain only two reachable rows');

end


