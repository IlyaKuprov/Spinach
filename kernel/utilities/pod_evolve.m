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
% The action callback prevents an additional automatic trace shift.
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

% Retain only coordinates on structural paths from the state to the detector
terms=model.terms; nranks=size(rho,1); adjacency=sparse(nranks,nranks);
for n=active, adjacency=adjacency+spones(sparse(terms{n})); end
coordinates=intersect(reachable_support(adjacency.',find(rho)),...
                      reachable_support(adjacency,find(coil)));
if isempty(coordinates), fid(:)=0; return; end
for n=active, terms{n}=terms{n}(coordinates,coordinates); end
[rho,rebase]=binary_normalise(rho(coordinates)); state_exp=state_exp+rebase;
[coil,rebase]=binary_normalise(coil(coordinates)); coil_exp=coil_exp+rebase;

% Determine a common binary frame for the complete time-scaled affine sum
[time_frac,time_exp]=log2(timestep); term_exp=zeros(size(coefficients));
coef_frac=complex(zeros(size(coefficients))); coef_exp=zeros(size(coefficients));
for n=active
    [coef_frac(n),coef_exp(n)]=binary_normalise(coefficients(n));
    largest=max(max(abs(real(terms{n})),[],'all'),max(abs(imag(terms{n})),[],'all'));
    if largest==0, term_exp(n)=-Inf; else, [~,term_exp(n)]=log2(largest); end
end
exponents=term_exp+coef_exp+time_exp; common_exp=max(exponents(active));
if ~isfinite(common_exp), fid(:)=fid(1); return; end

% Sum normalised terms before restoring their physical time-scaled magnitude
[~,order]=sort(exponents(active),'descend');
A=complex(sparse(numel(coordinates),numel(coordinates)));
for n=active(order)
    if ~isfinite(term_exp(n)), continue; end
    [term,~]=binary_normalise(terms{n});
    term=checked_product(term,checked_product(coef_frac(n),time_frac));
    scaled=complex(pow2(real(term),exponents(n)-common_exp),...
                   pow2(imag(term),exponents(n)-common_exp));
    if any((real(term)~=0)&(real(scaled)==0),'all')||...
       any((imag(term)~=0)&(imag(scaled)==0),'all')
        error('affine term dynamic range exceeds the normalised double representation.');
    end
    A=A+scaled;
end
scaled=complex(pow2(real(A),common_exp),pow2(imag(A),common_exp));
if any((real(A)~=0)&(real(scaled)==0),'all')||...
   any((imag(A)~=0)&(imag(scaled)==0),'all')||...
   any(isnan(scaled),'all')||any(isinf(scaled),'all')
    error('time-scaled generator must be representable in double precision.');
end
A=-1i*scaled; clear term scaled
if ~any(any(A,1)&(rho.'~=0))||~any(any(A,2)&(coil~=0)), fid(:)=fid(1); return; end

% Certify scalar physical underflow from the largest finite complex signal
if size(A,1)==1&&real(A)<log(realmin)+log(eps)-log(realmax)-2*log(2)
    fid(2:end)=0; return;
end

% Choose an equivalent action with the smaller initial coordinate growth rate
use_adjoint=real_growth(A',coil)<real_growth(A,rho);
if use_adjoint
    A=A'; [rho,coil]=deal(coil,rho); [state_exp,coil_exp]=deal(coil_exp,state_exp);
    fid(1)=conj(fid(1));
end

% Preserve the caller stream used by native action norm estimates
saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng)); completed=0;
while completed<nsteps

    % Separate the current state amplitude from a controlled real scalar shift
    growth=real_growth(A,rho); growth_exp=floor(growth/log(2));
    if ~isfinite(growth_exp)||abs(growth_exp)>=flintmax
        growth=0; growth_exp=0;
    end
    shifted=A-growth*speye(size(A,1)); width=min(size(A,1),nsteps-completed);

    % Bound action storage and growth without altering the physical generator
    hermitian=shifted/2+shifted'/2; diagonal=real(diag(hermitian));
    hermitian=hermitian-spdiags(diag(hermitian),0,size(A,1),size(A,1));
    growth_bound=max(diagonal+sum(abs(hermitian),2));
    if growth_bound>0
        width=min(width,max(1,floor((log(realmax)-log(norm(rho)))/growth_bound)));
    end

    % Evaluate only the requested actions with no implicit trace shift
    if size(A,1)==1
        actions=exp(shifted)*rho;
    else
        if ~isfinite(norm(shifted,1))||norm(shifted,1)>=flintmax
            error('working action norm exceeds the exact double step range.');
        end
        actions=expmv(@(flag,data)exponential_action(shifted,flag,data),rho,0:width);
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

% Estimate an amplitude shift without requiring a full physical-scale action
function growth=real_growth(A,rho)

% Any finite real scalar shift is valid; its estimate is not an invariance test
if size(A,1)==1
    growth=real(A);
else
    [scaled,exponent]=binary_normalise(A);
    growth=pow2(real(rho'*(scaled*rho))/(rho'*rho),exponent);
end

end

% Supply native exponential actions and adjoints without a second scalar shift
function result=exponential_action(A,flag,data)

% Follow the documented operator-action protocol used by expmv
switch flag
    case 'real', result=isreal(A);
    case 'notransp', result=checked_product(A,data);
    case 'transp', result=checked_product(A',data);
end

end

% Separate a finite real or complex array into bounded components and a binary scale
function [data,exponent]=binary_normalise(data)
largest=max(max(abs(real(data)),[],'all'),max(abs(imag(data)),[],'all'));
[~,exponent]=log2(largest);
scaled=complex(pow2(real(data),-exponent),pow2(imag(data),-exponent));
if any((real(data)~=0)&(real(scaled)==0),'all')||...
   any((imag(data)~=0)&(imag(scaled)==0),'all')||...
   any(isnan(scaled),'all')||any(isinf(scaled),'all')
    error('working-array dynamic range exceeds the normalised double representation.');
end
data=scaled;
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

% Find all directed paths from a support using one auxiliary graph vertex
function coordinates=reachable_support(adjacency,support)

% Traverse the structural graph without numerical actions or dense expansion
nranks=size(adjacency,1); connectivity=digraph(adjacency);
connectivity=addedge(connectivity,repmat(nranks+1,numel(support),1),support,ones(numel(support),1));
coordinates=bfsearch(connectivity,nranks+1); coordinates(coordinates==nranks+1)=[];

end

% Reject unsupported operand ranges before a bounded matrix product
function product=checked_product(left,right)

% Bound every nonzero real component product without forming a product array
left_parts=[nonzeros(real(left));nonzeros(imag(left))];
right_parts=[nonzeros(real(right));nonzeros(imag(right))];
if ~isempty(left_parts)&&~isempty(right_parts)&&min(abs(left_parts))*min(abs(right_parts))==0
    error('working product dynamic range exceeds the normalised double representation.');
end
product=left*right;

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


