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
% MATLAB expm supplies the maintained scaling/squaring exponential;
% see Al-Mohy and Higham, DOI:10.1137/09074721X. One propagator is
% amortised across the requested steps of the small reduced model.
% Coefficients, terms and time are assembled in one binary-exponent frame,
% restoring only the summed time-scaled generator. No fallback order is used.
% The propagator's real scalar growth and working-state amplitude are carried
% separately; only the detected observable is restored to physical units.
% This avoids amplitude overflow when a small detector gives a finite signal.
% Inputs and integer counts are double; nsteps must be below flintmax.
% Nonfinite projected data, unrepresentable generator/component ranges,
% and nonfinite output signals are explicit errors. No arbitrary-precision or
% componentwise accuracy guarantee is implied for ill-conditioned dynamics.
% Matrix products require nonzero operand-component pairs to remain nonzero
% at the working scale; this can conservatively reject a finite physical result.
% Exponent-aware component products preserve finite scaled observables.
% Structurally unreachable or unobservable coordinates are removed exactly.
% Constant requests and common structural stationary supports avoid expm;
% computed zero actions are never used to infer stationarity. Sparse terms
% remain sparse until the reduced dense exponential is actually needed.
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

% Choose a representable scalar growth factor before the maintained exponential
nranks=size(A,1); growth=sum(real(diag(A))/nranks);
growth_exp=floor(growth/log(2));
if ~isfinite(growth_exp)||abs(growth_exp)>=flintmax
    growth=0; growth_exp=0;
end
P=expm(full(A-growth*speye(nranks))); clear A
if any(isnan(P),'all')||any(isinf(P),'all')
    error('shifted reduced propagator must be finite.');
end
[P,prop_exp]=binary_normalise(P);
P=P*exp(growth-growth_exp*log(2)); prop_exp=prop_exp+growth_exp;

% Propagate normalised coordinates and restore only the observable amplitude
for n=1:nsteps
    [rho,step_exp]=binary_normalise(checked_product(P,rho)); state_exp=state_exp+prop_exp+step_exp;
    if abs(state_exp)>=flintmax, error('state exponent exceeds the exact integer range.'); end
    fid(n+1)=detect_scaled(coil,rho,state_exp+coil_exp);
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


