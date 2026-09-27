% Tests giant-spin cache identity in both insertion orders. Syntax:
%
%                    result=test_giant_cache()
%
% Outputs:
%
%     result - regression checks for retention, coefficients, and crystal
%
% Physical non-axial rank-two and rank-four terms have nonzero tensor and
% sample frames. Uncached production assembly supplies the references.
%
% talos@spindynamics.org

function result=test_giant_cache()

% Initialise the regression record
result=new_test_result('kernel/giant_cache','Giant-spin cache identity',...
                       'cached assembly must preserve giant coefficients and retention.');

% Specify physical Hermitian giant terms in rotated tensor frames
sys.magnet=0.73; sys.output='hush'; sys.disable={'hygiene'};
sys.enable={'ham_cache'}; sys.parallel={'processes',1}; sys.parprops={};
inter.giant.coeff={{[0 0 0],[2e5 3e5 7e5 -3e5 2e5],...
                    zeros(1,7),[4e4 0 2e4 0 8e4 0 2e4 0 4e4]}};
inter.giant.euler={{[0 0 0],[0.31 0.53 0.17],[0 0 0],[0.19 0.41 0.37]}};
bas.formalism='zeeman-hilb'; bas.approximation='none';
angles=[0.29 0.47 0.61]; isotopes={'E5','E6'};

% Exercise opposite cache insertion orders on distinct physical systems
for ordering=1:2

    % Build independently retained full and Zeeman-only systems
    sys.isotopes=isotopes(ordering);
    spin_system=basis(create(sys,inter),bas);
    full_system=assume(spin_system,'labframe');
    zeeman_system=assume(spin_system,'labframe','zeeman');
    systems={full_system,zeeman_system};
    references=cell(1,2); invariants=cell(1,2);

    % Obtain uncached references with identical physical inputs
    for n=1:2
        reference=systems{n}; reference.sys.enable={};
        [H,Q]=hamiltonian(reference);
        invariants{n}=H; references{n}=H+orientation(Q,angles);
    end

    % Require a Hermitian complex Hamiltonian with noncommuting giant terms
    assert(norm(references{1}-references{1}','fro')<1e-7);
    assert(norm(imag(references{1}),'fro')>1e5);
    assert(norm(references{1}-references{2},'fro')>1e6);
    assert(norm(references{1}*references{2}-references{2}*references{1},'fro')>1e12);

    % Compare both cache insertion orders and repeat cache hits
    for n=[ordering 3-ordering ordering]
        [H,Q]=hamiltonian(systems{n});
        result=test_close(result,sprintf('order %d retention %d',ordering,n),...
                          H+orientation(Q,angles),references{n},1e-7,1e-12,...
                          'giant retention must distinguish cache entries');
        result=test_close(result,'one-output invariant',hamiltonian(systems{n}),...
                          invariants{n},1e-7,1e-12,'output arity must preserve the invariant');
    end

    % Change only giant coefficients with the existing cache still populated
    changed=full_system;
    changed.inter.giant.coeff{1}{2}=1.3*changed.inter.giant.coeff{1}{2};
    reference=changed; reference.sys.enable={};
    [H,Q]=hamiltonian(reference); changed_ref=H+orientation(Q,angles);
    assert(norm(changed_ref-references{1},'fro')>1e5);
    [H,Q]=hamiltonian(changed);
    result=test_close(result,'changed coefficients',H+orientation(Q,angles),...
                      changed_ref,1e-7,1e-12,'giant coefficients must distinguish cache entries');

    % Confirm that ignored coefficients do not affect the Zeeman Hamiltonian
    changed=assume(changed,'labframe','zeeman');
    [H,Q]=hamiltonian(changed);
    result=test_close(result,'ignored coefficients',H+orientation(Q,angles),...
                      references{2},1e-7,1e-12,'Zeeman retention must ignore giant terms');

    % Exercise the production crystal full and Zeeman split
    parameters.spins=sys.isotopes; parameters.offset=0;
    parameters.orientation=angles; parameters.needs={'zeeman_op'};
    observed=crystal(spin_system,@cache_generators,parameters,'labframe');
    reference=spin_system; reference.sys.enable={};
    expected=crystal(reference,@cache_generators,parameters,'labframe');
    result=test_close(result,'crystal decomposition',observed,expected,1e-7,1e-12,...
                      'crystal must receive distinct full and Zeeman Hamiltonians');
end

% Preserve ordinary systems in both existing cache-key branches
for modes_present=0:1
    sys.isotopes={'1H'}; inter=struct(); inter.zeeman.scalar={2};
    if modes_present
        sys.isotopes={'1H','C3'}; inter.zeeman.scalar={2,[]};
        inter.modes.frqs={[],1000}; inter.modes.anharms={[],23};
    end
    spin_system=assume(basis(create(sys,inter),bas),'labframe');
    reference=spin_system; reference.sys.enable={};
    expected=hamiltonian(reference);
    for n=1:2
        result=test_close(result,'ordinary cache branch',hamiltonian(spin_system),...
                          expected,1e-7,1e-12,'non-giant cached assembly must be unchanged');
    end
end

end

% Capture both physical context generators for each insertion order
function answer=cache_generators(spin_system,parameters,H,R,K)

% Validate the context boundary arguments
grumble(spin_system,parameters,H,R,K);

% Return the full and Zeeman Hamiltonians without altering them
answer=[H parameters.hzeeman];

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system), error('spin_system must be a structure.'); end
if ~isstruct(parameters)||~isfield(parameters,'hzeeman')
    error('parameters must contain hzeeman.');
end
if ~isnumeric(H)||~ismatrix(H)||size(H,1)~=size(H,2)
    error('H must be a numeric square matrix.');
end
if ~isnumeric(R)||~isnumeric(K)||~isnumeric(parameters.hzeeman)||...
   ~isequal(size(H),size(R),size(K),size(parameters.hzeeman))
    error('all generators must have matching numeric dimensions.');
end
end


