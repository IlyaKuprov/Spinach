% Observable trajectory of a prepared fixed-basis POD model. Syntax:
%
%              fid=pod_evolve(model,coefficients,timestep,nsteps)
%
% Parameters:
%
%     model        - scalar structure from pod_project, containing a cell
%                    column of r-by-r Liouvillian terms, and r-by-1 rho/coil
%
%     coefficients - finite double column matching model.terms; complex
%                    coefficients are permitted, without inferred units
%
%     timestep     - finite nonnegative scalar time step, seconds
%
%     nsteps       - nonnegative integer number of evolution steps
%
% Outputs:
%
%     fid - complex 1-by-(nsteps+1) trace, including the initial observable
%
% Supplied terms use L=H+1i*R+1i*K in rad/s. Fixed-basis affine
% Galerkin evaluation follows Benner et al., DOI:10.1137/130932715.
% MATLAB exponential actions follow Al-Mohy and Higham,
% DOI:10.1137/100788860. Actions are streamed in blocks no wider than
% the reduced rank, never forming a full propagator or full state trajectory.
% Scalar models use the identical closed-form scalar exponential.
% Coefficients, terms and time are assembled in one binary-exponent frame,
% restoring only the summed time-scaled generator. No fallback order is used.
% A state-conditioned real scalar shift and working-state amplitude are carried
% separately; only the detected observable is restored to physical units.
% Equivalent forward/adjoint actions choose the smaller initial norm-growth
% rate; this changes neither the generator nor the Hermitian observable.
% The diagonal is stored separately; callbacks apply its scalar shift without
% constructing a shifted dense operator. The Hermitian bound and induced norms
% use sparse reductions or bounded dense scans, never a dense Hermitian part.
% Real actions and states retain real working storage. The native callback prevents an
% additional automatic trace shift.
% Block widths use a Hermitian-part growth bound and the double range;
% unscalarised action norms above flintmax are unsupported. Working actions
% that round entirely to zero are rejected, not declared stationary.
% Caller RNG is restored after native norm estimates, including on errors.
% A scalar damping bound certifies zero when every future signal component
% is below half the smallest subnormal, without testing a computed action.
% Physical signals below the subnormal range round to zero; nonfinite
% physical signals are errors.
% Inputs and integer counts are double; nsteps must be below flintmax.
% Nonfinite projected data, unrepresentable generator/component ranges,
% and nonfinite output signals are explicit errors. No arbitrary-precision or
% componentwise accuracy guarantee is implied for ill-conditioned dynamics.
% Matrix products require nonzero operand-component pairs to remain nonzero
% at the working scale; this can conservatively reject a finite physical result.
% Exponent-aware component products preserve finite scaled observables.
% Structurally unreachable or unobservable coordinates are removed exactly.
% Implicit support traversal uses bounded dense row/column slices; sparse
% reverse traversal stores transposed sparse terms. Fixed off-diagonal component
% bounds, Hermitian radii, and induced norms preserve the operator storage;
% shifted diagonal components are checked once per action block.
% Constant requests and common structural stationary supports avoid action evaluation;
% computed zero actions are never used to infer stationarity. Sparse terms
% remain sparse through the maintained exponential actions.
% Training, projection and evaluation remain separate. No physical system,
% metric, reference, processing policy or positivity constraint is inferred.
% Snapshot reconstruction energy does not certify observable accuracy.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pod_evolve.m>

function fid=pod_evolve(model,coefficients,timestep,nsteps)

% Validate the supplied prepared coordinates and evaluation grid
grumble(model,coefficients,timestep,nsteps);

% Separate coordinate amplitude from the working vectors before detection
[rho,state_exp]=binary_normalise(model.rho);
[coil,coil_exp]=binary_normalise(model.coil);
fid=complex(zeros(1,nsteps+1)); fid(1)=detect_scaled(coil,rho,state_exp+coil_exp);
if nsteps==0||timestep==0||~any(coefficients)||~any(rho)||~any(coil)
    fid(:)=fid(1); return;
end

% Certify common same-side structural stationarity before irrelevant scaling
active=find(coefficients).'; state_fixed=true; coil_fixed=true;
for n=active
    state_fixed=state_fixed&&~any(any(model.terms{n},1)&(rho.'~=0));
    coil_fixed=coil_fixed&&~any(any(model.terms{n},2)&(coil~=0));
end
if state_fixed||coil_fixed, fid(:)=fid(1); return; end

% Retain structural paths without materialising a dense support graph
terms=model.terms;
coordinates=intersect(reachable_support(terms,active,find(rho),false),...
                      reachable_support(terms,active,find(coil),true));
if isempty(coordinates), fid(:)=0; return; end
if numel(coordinates)<size(rho,1)
    for n=active, terms{n}=terms{n}(coordinates,coordinates); end
    [rho,rebase]=binary_normalise(rho(coordinates)); state_exp=state_exp+rebase;
    [coil,rebase]=binary_normalise(coil(coordinates)); coil_exp=coil_exp+rebase;
end

% Determine a common binary frame for the complete time-scaled affine sum
[time_frac,time_exp]=log2(timestep); term_exp=zeros(size(coefficients));
coef_frac=zeros(size(coefficients)); coef_exp=zeros(size(coefficients));
for n=active
    [coef_frac(n),coef_exp(n)]=binary_normalise(coefficients(n));
    largest=max(max(abs(real(terms{n})),[],'all'),max(abs(imag(terms{n})),[],'all'));
    if largest==0, term_exp(n)=-Inf; else, [~,term_exp(n)]=log2(largest); end
end
exponents=term_exp+coef_exp+time_exp; common_exp=max(exponents(active));
if ~isfinite(common_exp), fid(:)=fid(1); return; end

% Sum normalised terms before restoring their physical time-scaled magnitude
[~,order]=sort(exponents(active),'descend');
A=sparse(numel(coordinates),numel(coordinates));
for n=active(order)
    if ~isfinite(term_exp(n)), continue; end
    [term,~]=binary_normalise(terms{n});
    term=checked_product(term,checked_product(coef_frac(n),time_frac));
    scaled=binary_scale(term,exponents(n)-common_exp);
    A=A+scaled;
end

% Release affine term workspaces before restoring the physical generator
clear term scaled terms
scaled=binary_scale(A,common_exp);
A=-1i*scaled; clear scaled
if ~any(any(A,1)&(rho.'~=0))||~any(any(A,2)&(coil~=0)), fid(:)=fid(1); return; end

% Certify scalar physical underflow from the largest finite complex signal
if size(A,1)==1&&real(A)<log(realmin)+log(eps)-log(realmax)-2*log(2)
    fid(2:end)=0; return;
end

% Separate the diagonal once and cache bounded off-diagonal action estimates
diagonal=diag(A); A(1:size(A,1)+1:end)=0;
[radii,norms,minimum,scale_exp]=action_bounds(A,diagonal);

% Choose an equivalent action without transposing the stored generator
use_adjoint=real_growth(A,diagonal,scale_exp,coil)<real_growth(A,diagonal,scale_exp,rho);
if use_adjoint
    [rho,coil]=deal(coil,rho); [state_exp,coil_exp]=deal(coil_exp,state_exp);
    fid(1)=conj(fid(1));
end

% Preserve the caller stream used by native action norm estimates
saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng)); completed=0;
while completed<nsteps

    % Separate the current state amplitude from a controlled real scalar shift
    growth=real_growth(A,diagonal,scale_exp,rho); growth_exp=floor(growth/log(2));
    if ~isfinite(growth_exp)||abs(growth_exp)>=flintmax
        growth=0; growth_exp=0;
    end
    shifted_diag=diagonal-growth; width=min(size(A,1),nsteps-completed);

    % Bound action storage and growth without altering the physical generator
    growth_bound=max(real(shifted_diag)+radii);
    if growth_bound>0
        width=min(width,max(1,floor((log(realmax)-log(norm(rho)))/growth_bound)));
    end

    % Evaluate only the requested actions with no implicit trace shift
    if size(A,1)==1
        if use_adjoint, shifted_diag=conj(shifted_diag); end
        actions=exp(shifted_diag)*rho;
    else
        action_norm=max(norms(:,1+use_adjoint)+abs(shifted_diag));
        if ~isfinite(action_norm)||action_norm>=flintmax
            error('working action norm exceeds the exact double step range.');
        end
        component_bound=min(minimum,component_min(shifted_diag));
        actions=expmv(@(flag,data)exponential_action(A,shifted_diag,component_bound,...
                                                  use_adjoint,flag,data),rho,0:width);
        actions=actions(:,2:end);
    end

    % Restore each observable while carrying the last normalised state forward
    base_exp=state_exp;
    for k=1:width
        [column,column_exp]=binary_normalise(actions(:,k));
        if ~any(column), error('working exponential action underflowed to zero.'); end
        fractional=k*(growth-growth_exp*log(2)); fractional_exp=floor(fractional/log(2));
        column=column*exp(fractional-fractional_exp*log(2));
        [column,rebase]=binary_normalise(column);
        state_exp=base_exp+column_exp+k*growth_exp+fractional_exp+rebase;
        if abs(state_exp)>=flintmax, error('state exponent exceeds the exact integer range.'); end
        fid(completed+k+1)=detect_scaled(coil,column,state_exp+coil_exp);
    end
    rho=column; completed=completed+width;

end
if use_adjoint, fid=conj(fid); end
clear restore_rng

end

% Cache action bounds with sparse reductions or bounded dense workspace
function [radii,norms,minimum,exponent]=action_bounds(A,diagonal)

% Measure the Hermitian row radii and both induced column norms once
largest=max(max(abs(real(diagonal)),[],'all'),max(abs(imag(diagonal)),[],'all'));
if issparse(A)
    radii=full(sum(abs(A/2+A'/2),2));
    norms=full([sum(abs(A),1).' sum(abs(A),2)]); minimum=component_min(A);
    largest=full(max(largest,max(max(abs(real(A)),[],'all'),max(abs(imag(A)),[],'all'))));
else
    radii=zeros(size(A,1),1); norms=zeros(size(A,1),2); minimum=Inf;
    for k=1:size(A,1)
        column=A(:,k); row=A(k,:)';
        radii(k)=sum(abs(column/2+row/2)); norms(k,:)=[sum(abs(column)) sum(abs(row))];
        minimum=min(minimum,component_min(column));
        largest=max(largest,max(max(abs(real(column)),[],'all'),max(abs(imag(column)),[],'all')));
    end
end
[~,exponent]=log2(largest);

end

% Estimate the real Rayleigh rate without a full dense normalised generator
function growth=real_growth(A,diagonal,exponent,rho)

% The real Rayleigh rate is identical for a generator and its Hermitian adjoint
if issparse(A)
    action=binary_scale(A,-exponent)*rho;
else
    action=zeros(size(rho));
    for k=1:size(A,2)
        action=action+binary_scale(A(:,k),-exponent)*rho(k);
    end
end
action=action+binary_scale(diagonal,-exponent).*rho;
growth=pow2(real(rho'*action)/(rho'*rho),exponent);

end

% Apply shifted diagonal actions without copying or transposing the generator
function result=exponential_action(A,diagonal,minimum,adjoint,flag,data)

% Check changing operands and select the equivalent action or its adjoint
if strcmp(flag,'real'), result=isreal(A)&&isreal(diagonal); return; end
if minimum*component_min(data)==0
    error('working product dynamic range exceeds the normalised double representation.');
end
if xor(strcmp(flag,'transp'),adjoint)
    result=A'*data+conj(diagonal).*data;
else
    result=A*data+diagonal.*data;
end

end

% Separate a finite real or complex array into bounded components and a binary scale
function [data,exponent]=binary_normalise(data)
largest=max(max(abs(real(data)),[],'all'),max(abs(imag(data)),[],'all'));
[~,exponent]=log2(largest);
data=binary_scale(data,-exponent);
end

% Scale real or complex components without forcing complex working storage
function scaled=binary_scale(data,exponent)
scaled=pow2(real(data),exponent);
if ~isreal(data), scaled=complex(scaled,pow2(imag(data),exponent)); end
if nnz(real(data))~=nnz(real(scaled))||...
   (~isreal(data)&&nnz(imag(data))~=nnz(imag(scaled)))||...
   any(isnan(scaled),'all')||any(isinf(scaled),'all')
    error('working-array dynamic range exceeds the normalised double representation.');
end
end

% Restore an observable using exponent-aware real component products
function signal=detect_scaled(coil,rho,exponent)

% Accumulate the real and imaginary Hermitian products in independent frames
components=cell(2,2);
components(1,:)={real(coil),real(rho)};
components(2,:)={imag(coil),imag(rho)};
signal=complex(0);
for n=1:2
    fractions=[]; powers=[];
    for k=1:2
        left=components{k,1}; right=components{k,2}; indices=find((left~=0)&(right~=0));
        [left_frac,left_exp]=log2(abs(full(left(indices))));
        [right_frac,right_exp]=log2(abs(full(right(indices))));
        fractions=[fractions;sign(left(indices)).*sign(right(indices)).*left_frac.*right_frac]; %#ok<AGROW>
        powers=[powers;left_exp+right_exp]; %#ok<AGROW>
    end
    if ~isempty(powers)
        [powers,~,groups]=unique(powers); fractions=accumarray(groups,full(fractions));
        powers=powers(fractions~=0); fractions=fractions(fractions~=0);
        if ~isempty(powers)
            common=max(powers); scaled=pow2(fractions,powers-common);
            if any((fractions~=0)&(scaled==0))
                error('observable product dynamic range exceeds the normalised double representation.');
            end
            signal=signal+(1i^(n-1))*pow2(sum(scaled),common+exponent);
        end
    end
    components(1,:)={real(coil),imag(rho)};
    components(2,:)={-imag(coil),real(rho)};
end
if ~isfinite(signal), error('detected observable must be finite.'); end

end

% Traverse implicit term supports with bounded dense neighbour storage
function coordinates=reachable_support(terms,active,support,backward)

% Track each coordinate once without constructing an aggregate edge set
nranks=size(terms{active(1)},1); visited=false(nranks,1); queue=zeros(nranks,1);
tail=numel(support); queue(1:tail)=support; visited(support)=true; head=1;
if tail==nranks, coordinates=(1:nranks).'; return; end

% Retain column traversal for sparse adjoints without transposing dense terms
if backward
    for n=active
        if issparse(terms{n}), terms{n}=terms{n}.'; end
    end
end

% Discover neighbours from bounded dense slices or stored sparse columns
while head<=tail
    coordinate=queue(head); head=head+1;
    for n=active
        if backward&&~issparse(terms{n})
            neighbours=find(terms{n}(coordinate,:).');
        else
            neighbours=find(terms{n}(:,coordinate));
        end
        neighbours=neighbours(~visited(neighbours)); visited(neighbours)=true;
        queue(tail+(1:numel(neighbours)))=neighbours; tail=tail+numel(neighbours);
        if tail==nranks, coordinates=(1:nranks).'; return; end
    end
end
coordinates=find(visited);

end

% Reject unsupported operand ranges before a bounded matrix product
function product=checked_product(left,right)

% Bound every nonzero real component product without forming a product array
if component_min(left)*component_min(right)==0
    error('working product dynamic range exceeds the normalised double representation.');
end
product=left*right;

end

% Find the smallest nonzero component with bounded dense column storage
function minimum=component_min(data)

% Preserve stored sparse traversal without copying dense matrix components
minimum=Inf;
if issparse(data)
    parts=nonzeros(real(data));
    if ~isempty(parts), minimum=min(abs(parts)); end
    parts=nonzeros(imag(data));
    if ~isempty(parts), minimum=min(minimum,min(abs(parts))); end
else
    for k=1:size(data,2)
        parts=real(data(:,k)); parts=parts(parts~=0);
        if ~isempty(parts), minimum=min(minimum,min(abs(parts))); end
        parts=imag(data(:,k)); parts=parts(parts~=0);
        if ~isempty(parts), minimum=min(minimum,min(abs(parts))); end
    end
end

end

% Consistency enforcement
function grumble(model,coefficients,timestep,nsteps)
if ~isstruct(model)||~isscalar(model)||~all(isfield(model,{'terms','rho','coil'}))
    error('model must contain terms, rho, and coil.');
end
if ~iscell(model.terms)||~iscolumn(model.terms)||isempty(model.terms)
    error('model.terms must be a nonempty cell column.');
end
nranks=size(model.terms{1},1);
for n=1:numel(model.terms)
    if ~isa(model.terms{n},'double')||nranks<1||...
       ~isequal(size(model.terms{n}),[nranks nranks])||...
       any(isnan(model.terms{n}),'all')||any(isinf(model.terms{n}),'all')
        error('model terms must be conformable finite square double matrices.');
    end
end
if ~isa(model.rho,'double')||~isequal(size(model.rho),[nranks 1])||...
   any(isnan(model.rho),'all')||any(isinf(model.rho),'all')||~isa(model.coil,'double')||...
   ~isequal(size(model.coil),[nranks 1])||any(isnan(model.coil),'all')||any(isinf(model.coil),'all')
    error('model rho and coil must be conformable finite double columns.');
end
if ~isa(coefficients,'double')||~isequal(size(coefficients),size(model.terms))||...
   ~all(isfinite(coefficients))
    error('coefficients must be a matching finite double column.');
end
if ~isa(timestep,'double')||~isreal(timestep)||~isscalar(timestep)||...
   ~isfinite(timestep)||timestep<0
    error('timestep must be a finite nonnegative scalar.');
end
if ~isa(nsteps,'double')||~isreal(nsteps)||~isscalar(nsteps)||...
   ~isfinite(nsteps)||nsteps<0||nsteps>=flintmax||mod(nsteps,1)~=0
    error('nsteps must be a nonnegative integer.');
end
end


