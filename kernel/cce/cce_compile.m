% Compile canonical tensors and metadata for native cluster expansion.
% Syntax:
%
%           model=cce_compile(spin_system,parameters,assumptions)
%
% spin_system is create output; parameters obeys cce_options. Tensors
% remain in rad/s and coordinates in Angstrom. No full-bath basis is built.
% The calling context resolves cce.numerics.atol, rtol, and max_memory_bytes.
% Observed and bath IDs retain the caller's order; groups retain physical IDs.
% A supplied basis selects zeeman-hilb or zeeman-liouv; absent basis selects
% zeeman-hilb. Stale full-system basis arrays are removed. Local representations
% are complete, not restricted bases. The model owns canonical tensors,
% assumptions, native parameters, isotope labels, dimensions, and projectors.
%
% Parameters:
%
%   spin_system - physical system returned by create
%   parameters  - native CCE options with resolved numerical tolerances,
%                 memory budget, and orientation in radians
%   assumptions - stock interaction assumption string or argument cell
%
% Outputs:
%
%   model - canonical angular-frequency tensors, physical-spin ordering,
%           coordinates in Angstrom, and retained approximation metadata
%
% ilya.kuprov@weizmann.ac.il

function model=cce_compile(spin_system,parameters,assumptions)

% Resolve native defaults without translating obsolete parameter forms
[parameters,~]=cce_options(spin_system,parameters,@cce_compile,assumptions);

% Validate the compiler inputs before copying or rotating physical tensors
grumble(spin_system,parameters,assumptions);

% Carry only a supported representation, never full-system basis arrays
formalism='zeeman-hilb';
if isfield(spin_system,'bas')
    formalism=spin_system.bas.formalism;
    spin_system=rmfield(spin_system,'bas');
end

% Apply the requested stock approximation without building a basis
if ischar(assumptions)
    assumptions={assumptions};
end
spin_system=assume(spin_system,assumptions{:});

% Copy immutable source metadata and stable physical identifiers
model.spin_system=spin_system;
model.parameters=parameters;
model.kind='cce_model';
model.formalism=formalism;
model.isotopes=spin_system.comp.isotopes;
model.groups={};
if isfield(parameters.cce,'groups')
    model.groups=parameters.cce.groups;
end
model.manifold_projectors=struct();
model.physical_ids=1:spin_system.comp.nspins;
model.observed_ids=parameters.observed;
model.bath_ids=parameters.cce.bath_spins;
model.multiplicities=spin_system.comp.mults;
model.gamma=spin_system.inter.gammas;
model.units=struct('frequency','rad/s','distance','Angstrom','angle','radian');
model.numerical_cutoffs=struct('interaction_Hz',spin_system.tols.inter_cutoff,...
                              'operator_rad_s',spin_system.tols.liouv_zero);
model.assumptions=assumptions;
model.orientation=parameters.orientation;
model.interaction_cutoff=inf;
if isfield(parameters.cce,'interaction_cutoff')
    model.interaction_cutoff=parameters.cce.interaction_cutoff;
end
model.numerics=parameters.cce.numerics;
model.backend=parameters.cce.backend;

% Preserve absent coordinates rather than inventing spatial positions
nspins=spin_system.comp.nspins;
model.coordinates_A=nan(nspins,3);
if isfield(spin_system.inter,'coordinates')
    for n=1:nspins
        if ~isempty(spin_system.inter.coordinates{n})
            model.coordinates_A(n,:)=spin_system.inter.coordinates{n};
        end
    end
end

% Rotate canonical Cartesian tensors into the requested laboratory frame
rotation=euler2dcm(parameters.orientation);
model.zeeman=cell(1,nspins);
model.quadratic=repmat({zeros(3)},1,nspins);
model.quadratic_strength=repmat({'ignore'},1,nspins);
model.pairs=struct('ids',{},'tensor',{},'retained',{});
model.omitted_pairs=zeros(0,2);
for n=1:nspins
    model.zeeman{n}=rotation*spin_system.inter.zeeman.matrix{n}*rotation.';
end

% Index the directed records and rotate the on-site terms
list=spin_system.inter.coupling.matrix;
from=cell2mat(list(:,1));
to=cell2mat(list(:,2));
index=sparse(from,to,1:size(list,1),nspins,nspins);
omit=false(size(list,1),1);
for r=find(from==to).'
    n=from(r);
    tensor=list{r,3};
    if norm(tensor,2)>2*pi*spin_system.tols.inter_cutoff
        model.quadratic{n}=rotation*tensor*rotation.';
    end
    model.quadratic_strength{n}=spin_system.inter.coupling.strength{r,3};
end

% Combine opposite directions in ascending unordered physical pair order
[k_ids,n_ids]=find(tril(index+index.',-1));
for p=1:numel(n_ids)

    % Accumulate both independently labelled directions
    n=n_ids(p);
    k=k_ids(p);
    tensor=zeros(3);
    retained=zeros(3);
    for direction=1:2

        % Retrieve the directed coupling and its aligned label
        if direction==1
            row=index(n,k);
        else
            row=index(k,n);
        end
        if row==0
            continue
        end
        value=list{row,3};
        if norm(value,2)<=2*pi*spin_system.tols.inter_cutoff
            continue
        end

        % Separate the stock isotropic and anisotropic numerical cutoffs
        isotropic=trace(value)/3;
        anisotropic=value-isotropic*eye(3);
        if abs(isotropic)<=2*pi*spin_system.tols.inter_cutoff
            isotropic=0;
        end
        [~,rank_one,rank_two]=mat2sphten(anisotropic);
        if norm(rank_one,2)<=spin_system.tols.liouv_zero&&...
           norm(rank_two,2)<=spin_system.tols.liouv_zero
            anisotropic=zeros(3);
        end
        value=isotropic*eye(3)+anisotropic;

        % Stock spherical pair convention transposes each ordered cell
        value=rotation*value.'*rotation.';
        part=retain_pair(value,spin_system.inter.coupling.strength{row,3});
        if direction==2
            value=value.';
            part=part.';
        end
        tensor=tensor+value;
        retained=retained+part;

    end

    % Apply only the explicitly requested physical distance cutoff
    if isfinite(model.interaction_cutoff)&&(norm(tensor,'fro')>0)
        if any(~isfinite(model.coordinates_A([n k],:)),'all')
            error('finite interaction_cutoff requires coordinates for every coupled pair.');
        end
        if norm(model.coordinates_A(n,:)-model.coordinates_A(k,:))>model.interaction_cutoff
            model.omitted_pairs(end+1,:)=[n k];
            rows=[index(n,k) index(k,n)];
            omit(rows(rows>0))=true;
            tensor=zeros(3);
            retained=zeros(3);
        end
    end

    % A cancelled bare tensor can still have a nonzero retained part
    if any(tensor,'all')||any(retained,'all')
        model.pairs(end+1)=struct('ids',[n k],'tensor',tensor,'retained',retained);
    end

end

% Remove pruned directions and their aligned strength records
spin_system.inter.coupling.matrix(omit,:)=[];
spin_system.inter.coupling.strength(omit,:)=[];

% Bind the physically pruned stock fallback to the same model
model.spin_system=spin_system;
model.canonical_roundoff=0;
for n=1:nspins

    % Bound tensor rotation and carrier subtraction before small residuals
    scale=norm(spin_system.inter.zeeman.matrix{n},'fro')+...
          abs(spin_system.inter.basefrqs(n));
    model.canonical_roundoff=model.canonical_roundoff+27*eps*scale;

end

% Exclude execution-only options from the physical model fingerprint
physical=parameters;
for name={'parallel','checkpoint','retain_raw','convergence'}
    if isfield(physical.cce,name{1})
        physical.cce=rmfield(physical.cce,name{1});
    end
end
model.source_hash=md5_hash({spin_system.comp,spin_system.inter,spin_system.tols,...
                            physical,assumptions,formalism});

% Declare the validated direct domain rather than guessing new physics
model.direct_supported=ismember(assumptions{1},{'labframe','nmr','esr','deer','deer-zz'});
model.fallback_reason='';
if ~model.direct_supported
    model.fallback_reason='unvalidated assumption set';
end
if any(~cellfun(@isempty,spin_system.inter.giant.coeff))
    model.direct_supported=false;
    model.fallback_reason='giant-spin spherical terms';
end
if isfield(spin_system.inter,'modes')
    model.direct_supported=false;
    model.fallback_reason='spin-boson terms';
end
if strcmp(model.backend,'direct')&&~model.direct_supported
    error('direct Hamiltonian is not validated for %s.',model.fallback_reason);
end

% Resolve the observed transition after compiling its native tensor model
if isfield(parameters,'transition')
    [model.transition,model.transition_meta]=cce_transition(model,parameters,assumptions);
    dimension=numel(model.transition.alpha);
    if 32*dimension^2>parameters.cce.numerics.max_memory_bytes
        error('observed projectors exceed numerics.max_memory_bytes.');
    end
    model.manifold_projectors.alpha=model.transition.alpha*model.transition.alpha';
    model.manifold_projectors.beta=model.transition.beta*model.transition.beta';
end

end

% Apply the stock Cartesian retention convention to an ordered pair
function tensor=retain_pair(tensor,strength)
switch strength
    case 'strong'

        % Retain every Cartesian component
    case 'secular'

        % Retain zero total Zeeman quantum order including antisymmetric xy
        value=(tensor(1,1)+tensor(2,2))/2;
        cross=(tensor(1,2)-tensor(2,1))/2;
        tensor=[value cross 0;-cross value 0;0 0 tensor(3,3)];
    case 'zz'

        % Retain the longitudinal product only
        tensor=diag([0 0 tensor(3,3)]);
    case 'z*'

        % Retain the first spin longitudinal components
        tensor(1:2,:)=0;
    case '*z'

        % Retain the second spin longitudinal components
        tensor(:,1:2)=0;
    case 'ignore'

        % Remove the entire ordered interaction
        tensor=zeros(3);
    otherwise

        % Leave unsupported domains for the stock backend
        tensor=zeros(3);

end

end

% Validate the compiler-owned geometry, representation, and numerical policy
function grumble(spin_system,parameters,assumptions)
if ~isstruct(spin_system)||~isscalar(spin_system)||~isfield(spin_system,'comp')||...
   ~isstruct(parameters)||~isscalar(parameters)||...
   ~all(isfield(parameters,{'observed','orientation','cce'}))
    error('cce_compile:input','Supply a Spinach system and native CCE parameters.');
end
if ~(ischar(assumptions)&&isrow(assumptions)&&~isempty(assumptions))&&...
   ~(iscell(assumptions)&&isrow(assumptions)&&~isempty(assumptions)&&ischar(assumptions{1}))
    error('cce_compile:assumptions','Supply explicit stock assume arguments.');
end
if ~isnumeric(parameters.orientation)||~isreal(parameters.orientation)||...
   ~isequal(size(parameters.orientation),[1 3])||any(~isfinite(parameters.orientation))
    error('cce_compile:orientation','Supply a finite 1-by-3 orientation in radians.');
end
for ids={parameters.observed,parameters.cce.bath_spins}
    value=ids{1};
    if ~isnumeric(value)||~isreal(value)||(~isempty(value)&&~isrow(value))||...
       any(~isfinite(value)|value<1|value>spin_system.comp.nspins|value~=fix(value))||...
       numel(unique(value))~=numel(value)
        error('cce_compile:spins','Supply distinct physical spin identifiers in rows.');
    end
end
if isempty(parameters.observed)||~isempty(intersect(parameters.observed,parameters.cce.bath_spins))
    error('cce_compile:spins','Observed spins must be nonempty and disjoint from the bath.');
end
if ~all(isfield(parameters.cce.numerics,{'atol','rtol','max_memory_bytes'}))
    error('cce_compile:policy','Resolve atol, rtol and max_memory_bytes before compilation.');
end
for name={'atol','rtol','max_memory_bytes'}
    value=parameters.cce.numerics.(name{1});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<0
        error('cce_compile:policy','Numerical budgets must be finite nonnegative scalars.');
    end
end
if parameters.cce.numerics.max_memory_bytes==0
    error('cce_compile:policy','The memory budget must be positive.');
end
if isfield(spin_system,'bas')&&...
   ~ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
    error('cce_compile:formalism','Only complete Zeeman representations are supported.');
end
end


