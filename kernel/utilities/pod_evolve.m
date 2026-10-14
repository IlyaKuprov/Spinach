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
% The assembled L is in rad/s: L=H+1i*R+1i*K. This function uses MATLAB
% expm(-1i*L*timestep), then repeats the fixed propagator. Training and
% projection are separate costs. No physical spin_system is modified, and
% no tolerance-based screening, processing, or reference is inferred.
% Accuracy is limited by the supplied reduction and generator, not by a
% snapshot-energy guarantee. No positivity or trace preservation is implied.
% Zero-step, zero-duration, and all-zero-coefficient requests validate inputs
% and return the initial signal, repeated as needed, without assembling a
% generator. A zero initial state or detector also needs no propagation.
% Zero coefficients are omitted from assembly; sparse active terms retain a
% sparse generator until a dense exponential is needed. An exactly zero
% assembled generator returns the constant initial signal without an exponential.
% Zero generator columns on the state support, or zero rows on the detector
% support, also give constant traces. Computed zero actions are not inferred.
% Sparse model terms and states are validated without expanding implicit zeros.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pod_evolve.m>

function fid=pod_evolve(model,coefficients,timestep,nsteps)

% Check the prepared model and evaluation contract
grumble(model,coefficients,timestep,nsteps);

% Return known constant signals without constructing a propagator
rho=model.rho; fid=complex(zeros(1,nsteps+1)); fid(1)=model.coil'*rho;
if (nsteps==0)||(timestep==0)||~any(coefficients)||~any(rho)||~any(model.coil)
    fid(:)=fid(1); return;
end

% Assemble the projected affine Liouvillian
L=complex(sparse(size(model.terms{1},1),size(model.terms{1},2)));
for n=find(coefficients).'
    L=L+coefficients(n)*model.terms{n};
end
if any(isnan(L),'all')||any(isinf(L),'all')
    error('assembled Liouvillian must be finite.');
end

% Return the constant observable for a structurally stationary state or detector
if ~any(any(L,1)&(rho.'~=0))||~any(any(L,2)&(model.coil~=0)), fid(:)=fid(1); return; end

% Propagate the reduced observable with the Spinach sign convention
P=expm(-1i*full(L)*timestep);
for n=1:nsteps
    rho=P*rho; fid(n+1)=model.coil'*rho;
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
if ~isnumeric(timestep)||~isreal(timestep)||~isscalar(timestep)||...
   ~isfinite(timestep)||timestep<0
    error('timestep must be a finite nonnegative scalar.');
end
if ~isnumeric(nsteps)||~isreal(nsteps)||~isscalar(nsteps)||...
   ~isfinite(nsteps)||nsteps<0||mod(nsteps,1)~=0
    error('nsteps must be a nonnegative integer.');
end
end


