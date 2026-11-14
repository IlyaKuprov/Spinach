% Compiles canonical Spinach metadata for cluster Hamiltonians. Syntax:
%
%              model=cce_prepare(spin_system,parameters)
%
% spin_system is returned by create. parameters supplies central_spins,
% bath_spins, assumptions (name or {name,retention}), orientation (rad),
% interaction_cutoff (Angstrom), and numerics.atol/rtol/max_memory_bytes.
% Optional hamiltonian_backend is auto, direct, or stock (default auto).
% The model stores rad/s tensors; isotope conversion is never repeated.
%
% Inputs:
%   spin_system - output of create, not unconverted sys/inter.
%   parameters - central_spins, bath_spins, assumptions (name or
%   {name,retention}), Euler orientation in rad, interaction_cutoff in
%   Angstrom, and numerics.atol, rtol, max_memory_bytes; optional central and
%   hamiltonian_backend.
%
% Outputs:
%   model - canonical source metadata, physical IDs and maps, local and pair
%   interactions, central projectors when supplied, backend capability and
%   content hash.
%
% See also cce_hamiltonian, cce_project.
%
% ilya.kuprov@weizmann.ac.il

function model=cce_prepare(spin_system,parameters)

% Check canonical metadata and explicit physical choices
grumble(spin_system,parameters);

% Apply the requested stock approximation without building a basis
assumptions=parameters.assumptions;
if ischar(assumptions), assumptions={assumptions}; end
spin_system=assume(spin_system,assumptions{:});

% Copy immutable source metadata and stable physical identifiers
model.spin_system=spin_system; model.parameters=parameters;
model.central_projectors=struct();
if isfield(parameters,'central')&&all(isfield(parameters.central,{'alpha','beta'}))
    dimension=numel(parameters.central.alpha);
    if 32*dimension^2>parameters.numerics.max_memory_bytes
        error('central projectors exceed numerics.max_memory_bytes.');
    end
    model.central_projectors.alpha=parameters.central.alpha*parameters.central.alpha';
    model.central_projectors.beta=parameters.central.beta*parameters.central.beta';
end
model.physical_ids=1:spin_system.comp.nspins;
model.central_ids=parameters.central_spins;
model.bath_ids=parameters.bath_spins;
model.multiplicities=spin_system.comp.mults;
model.gamma=spin_system.inter.gammas;
model.units=struct('frequency','rad/s','distance','Angstrom','angle','radian');
model.numerical_cutoffs=struct('interaction_Hz',spin_system.tols.inter_cutoff,...
                              'operator_rad_s',spin_system.tols.liouv_zero);
model.assumptions=assumptions; model.orientation=parameters.orientation;
model.interaction_cutoff=parameters.interaction_cutoff;
model.numerics=parameters.numerics; model.backend='auto';
if isfield(parameters,'hamiltonian_backend')
    model.backend=parameters.hamiltonian_backend;
end

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
model.zeeman=cell(1,nspins); model.quadratic=cell(1,nspins);
model.pairs=struct('ids',{},'tensor',{});
model.retained=cell(nspins); model.omitted_pairs=zeros(0,2);
for n=1:nspins

    % Rotate the Zeeman and on-site quadratic tensors
    model.zeeman{n}=rotation*spin_system.inter.zeeman.matrix{n}*rotation.';
    tensor=spin_system.inter.coupling.matrix{n,n};
    if isempty(tensor), tensor=zeros(3); end
    if norm(tensor,2)<=2*pi*spin_system.tols.inter_cutoff, tensor=zeros(3); end
    model.quadratic{n}=rotation*tensor*rotation.';

    % Combine ordered cells only after applying their retention masks
    for k=n+1:nspins
        tensor=zeros(3); retained=zeros(3);
        for direction=1:2

            % Retrieve one ordered cell and its independently declared strength
            ids=[n k]; if direction==2, ids=fliplr(ids); end
            value=spin_system.inter.coupling.matrix{ids(1),ids(2)};
            if isempty(value)||norm(value,2)<=2*pi*spin_system.tols.inter_cutoff
                continue
            end

            % Separate the stock isotropic and anisotropic numerical cutoffs
            isotropic=trace(value)/3; anisotropic=value-isotropic*eye(3);
            if abs(isotropic)<=2*pi*spin_system.tols.inter_cutoff, isotropic=0; end
            [~,rank_one,rank_two]=mat2sphten(anisotropic);
            if norm(rank_one,2)<=spin_system.tols.liouv_zero&&...
               norm(rank_two,2)<=spin_system.tols.liouv_zero
                anisotropic=zeros(3);
            end
            value=isotropic*eye(3)+anisotropic;

            % Stock spherical pair convention transposes each ordered cell
            value=rotation*value.'*rotation.';
            part=retain_pair(value,spin_system.inter.coupling.strength{ids(1),ids(2)});
            if direction==2, value=value.'; part=part.'; end
            tensor=tensor+value; retained=retained+part;
        end

        % Apply only the explicitly requested physical distance cutoff
        if isfinite(model.interaction_cutoff)&&(norm(tensor,'fro')>0)
            if any(~isfinite(model.coordinates_A([n k],:)),'all')
                error('finite interaction_cutoff requires coordinates for every coupled pair.');
            end
            if norm(model.coordinates_A(n,:)-model.coordinates_A(k,:))>model.interaction_cutoff
                model.omitted_pairs(end+1,:)=[n k];
                spin_system.inter.coupling.matrix{n,k}=zeros(3);
                spin_system.inter.coupling.matrix{k,n}=zeros(3);
                tensor=zeros(3); retained=zeros(3);
            end
        end

        % Keep a compact unordered pair list and a symmetric lookup
        if norm(tensor,'fro')>0
            model.pairs(end+1)=struct('ids',[n k],'tensor',tensor);
        end
        model.retained{n,k}=retained; model.retained{k,n}=retained.';
    end
end

% Bind the physically pruned stock fallback to the same model
model.spin_system=spin_system;
model.canonical_roundoff=0;
for n=1:nspins

    % Bound tensor rotation and carrier subtraction before small residuals
    scale=norm(spin_system.inter.zeeman.matrix{n},'fro')+...
          abs(spin_system.inter.basefrqs(n));
    model.canonical_roundoff=model.canonical_roundoff+27*eps*scale;
end
model.source_hash=md5_hash({spin_system.comp,spin_system.inter,parameters});

% Declare the validated direct domain rather than guessing new physics
model.direct_supported=ismember(assumptions{1},{'labframe','nmr','esr','deer','deer-zz'});
model.fallback_reason='';
if ~model.direct_supported, model.fallback_reason='unvalidated assumption set'; end
if any(~cellfun(@isempty,spin_system.inter.giant.coeff))
    model.direct_supported=false; model.fallback_reason='giant-spin spherical terms';
end
if isfield(spin_system.inter,'modes')
    model.direct_supported=false; model.fallback_reason='spin-boson terms';
end
if strcmp(model.backend,'direct')&&~model.direct_supported
    error('direct Hamiltonian is not validated for %s.',model.fallback_reason);
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

% Validate public inputs before creating any cluster-sized objects
function grumble(spin_system,parameters)
if ~isstruct(spin_system)||~isfield(spin_system,'comp')||...
   ~isfield(spin_system,'inter')||~isfield(spin_system.inter,'gammas')
    error('spin_system must be canonical metadata returned by create.');
end
if ~isstruct(parameters)||~isscalar(parameters)
    error('parameters must be a scalar structure.');
end
known={'central_spins','bath_spins','method','expansion','quantity','times',...
       'order','central','assumptions','orientation','hamiltonian_backend',...
       'interaction_cutoff','bath_radius','graph','groups','partition',...
       'bath_state','sampling','exterior','pulses','lindblad','second_order',...
       'numerics','convergence','output','execution','disjoint','interlaced',...
       'noise','observable','operators','controls','normalized','meanfield_radius',...
       'clusters','cluster_cap','meanfield','interlace','sequence','sequence_name',...
       'deer','jumps','reducer','parallel','workers','retain_raw','checkpoint',...
       'spectrum','analysis','averaging'};
if ~all(ismember(fieldnames(parameters),known))
    error('unknown CCE parameter field.');
end
required={'central_spins','bath_spins','assumptions','orientation',...
          'interaction_cutoff','numerics'};
if ~all(isfield(parameters,required)), error('missing mandatory model parameters.'); end
for name={'central_spins','bath_spins'}
    ids=parameters.(name{1});
    if ~isnumeric(ids)||~isreal(ids)||(~isempty(ids)&&~isrow(ids))||...
       any(~isfinite(ids))||any(ids~=fix(ids))||any(ids<1)||...
       any(ids>spin_system.comp.nspins)||numel(unique(ids))~=numel(ids)
        error('%s must contain unique physical spin indices.',name{1});
    end
end
if isempty(parameters.central_spins)||...
   ~isempty(intersect(parameters.central_spins,parameters.bath_spins))
    error('central_spins must be nonempty and disjoint from bath_spins.');
end
if isfield(parameters,'central')
    if ~isstruct(parameters.central)||~isscalar(parameters.central)
        error('central must be a scalar structure.');
    end
    if xor(isfield(parameters.central,'alpha'),isfield(parameters.central,'beta'))
        error('central.alpha and beta must be supplied together.');
    end
    if isfield(parameters.central,'alpha')
        dimension=prod(spin_system.comp.mults(parameters.central_spins));
        for name={'alpha','beta'}
            state=parameters.central.(name{1});
            if ~isa(state,'double')||~isequal(size(state),[dimension 1])||...
               any(~isfinite(state))||abs(norm(state)-1)>dimension*eps
                error('central vectors must be finite normalized double columns.');
            end
        end
        if abs(parameters.central.alpha'*parameters.central.beta)>dimension*eps
            error('central vectors must be orthogonal.');
        end
    end
end
assumptions=parameters.assumptions;
if ischar(assumptions), assumptions={assumptions}; end
if ~iscell(assumptions)||~ismember(numel(assumptions),[1 2])||...
   ~all(cellfun(@ischar,assumptions))
    error('assumptions must be a name or {name,retention}.');
end
if numel(assumptions)==2&&~ismember(assumptions{2},{'zeeman','couplings'})
    error('retention must be zeeman or couplings.');
end
if ~isnumeric(parameters.orientation)||~isreal(parameters.orientation)||...
   ~isequal(size(parameters.orientation),[1 3])||any(~isfinite(parameters.orientation))
    error('orientation must be a finite real 1-by-3 Euler vector in radians.');
end
cutoff=parameters.interaction_cutoff;
if ~isnumeric(cutoff)||~isreal(cutoff)||~isscalar(cutoff)||isnan(cutoff)||cutoff<0
    error('interaction_cutoff must be nonnegative, in Angstrom, or Inf.');
end
if isfield(parameters,'hamiltonian_backend')&&...
   ~ismember(parameters.hamiltonian_backend,{'auto','direct','stock'})
    error('hamiltonian_backend must be auto, direct, or stock.');
end
num=parameters.numerics;
if ~isstruct(num)||~all(isfield(num,{'atol','rtol','max_memory_bytes'}))
    error('numerics requires atol, rtol, and max_memory_bytes.');
end
for name={'atol','rtol','max_memory_bytes'}
    value=num.(name{1});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<0
        error('numerics.%s must be a finite nonnegative scalar.',name{1});
    end
end
if num.atol+num.rtol==0||num.max_memory_bytes<1||...
   num.max_memory_bytes~=fix(num.max_memory_bytes)
    error('nonzero error budget and positive integer memory ceiling required.');
end
if isfield(num,'precision')&&~strcmp(num.precision,'double')
    error('only double precision is supported.');
end
for n=1:spin_system.comp.nspins
    value=spin_system.inter.zeeman.matrix{n};
    if ~isequal(size(value),[3 3])||~isreal(value)||any(~isfinite(value),'all')
        error('canonical Zeeman tensors must be finite real 3-by-3 matrices.');
    end
    for k=1:spin_system.comp.nspins
        value=spin_system.inter.coupling.matrix{n,k};
        if ~isempty(value)&&(~isequal(size(value),[3 3])||...
           ~isreal(value)||any(~isfinite(value),'all'))
            error('canonical pair tensors must be finite real 3-by-3 matrices.');
        end
    end
    value=spin_system.inter.coupling.matrix{n,n};
    if ~isempty(value)&&norm(value,'fro')>0
        bound=9*eps*norm(value,'fro');
        if norm(value-value.','fro')>bound||abs(trace(value))>bound
            error('on-site quadratic tensors must be traceless and symmetric.');
        end
    end
end
end

