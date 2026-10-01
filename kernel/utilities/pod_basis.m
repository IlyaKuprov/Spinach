% Snapshot proper orthogonal decomposition in Euclidean state coordinates.
% Accepts snapshots or pre-truncated snapshot factors; supplied factors
% define the approximated data, not any omitted original trajectories.
% Syntax:
%
%              result=pod_basis(blocks,weights,target_rank,options)
%
% Parameters:
%
%     blocks      - nonempty cell column of finite double n-by-m matrices,
%                   with a common row count; complex entries are supported
%
%     weights     - matching real nonnegative double column; block k is
%                   multiplied by sqrt(weights(k)), without normalisation
%
%     target_rank - positive integer not exceeding the numerical data rank
%
%     options     - scalar structure with all four required fields:
%                   method: 'svd' or 'randomised'; oversampling: nonnegative
%                   integer; power_iters: nonnegative integer; seed: integer
%                   in [0,2^32-1]. For 'svd', oversampling and power_iters
%                   must be zero and seed is unused. For 'randomised', rank
%                   plus oversampling must not exceed either data dimension
%
% Outputs:
%
%     result.basis         - n-by-target_rank orthonormal complex basis
%     result.singular_vals - singular values of supplied weighted data for
%                            'svd', or of its range sketch for 'randomised'
%     result.total_energy  - squared Frobenius norm of weighted supplied data
%     result.residual_abs  - ||data-basis*(basis'*data)||_F, on supplied data
%     result.residual_rel  - residual_abs/sqrt(total_energy)
%     result.orth_error    - ||basis'*basis-eye(target_rank)||_F
%
% Exact SVD minimises the supplied weighted snapshot Frobenius projection
% error. Randomised sketches and pre-truncated factors qualify that claim.
% Residuals are measured data diagnostics, not observable error bounds.
% Rank-deficient inputs are rejected if the requested rank exceeds the
% numerical rank; no arbitrary null-space directions are supplied. The
% caller's random state is preserved, including on an error.
%
% ilya.kuprov@weizmann.ac.il

function result=pod_basis(blocks,weights,target_rank,options)

% Check the supplied data contract
grumble(blocks,weights,target_rank,options);

% Assemble explicitly weighted snapshot blocks
for n=1:numel(blocks)
    blocks{n}=sqrt(weights(n))*blocks{n};
end
data=full([blocks{:}]); total_energy=norm(data,'fro')^2;
if ~all(isfinite(data),'all')||~isfinite(total_energy)||total_energy==0
    error('weighted data must have finite nonzero energy.');
end

% Compute the exact basis or a seeded randomised range sketch
if strcmp(options.method,'svd')
    [modes,values,~]=svd(data,'econ');
else
    saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng));
    rng(options.seed,'twister'); sketch_rank=target_rank+options.oversampling;
    omega=complex(randn(size(data,2),sketch_rank),randn(size(data,2),sketch_rank));
    [range,~]=qr(data*omega,0);
    for n=1:options.power_iters
        [dual_range,~]=qr(data'*range,0);
        [range,~]=qr(data*dual_range,0);
    end
    [left_modes,values,~]=svd(range'*data,'econ'); modes=range*left_modes;
    clear restore_rng
end

% Reject unsupported rank rather than returning arbitrary null-space modes
singular_vals=diag(values);
rank_tol=max(size(data))*eps(singular_vals(1));
if singular_vals(target_rank)<=rank_tol
    error('target_rank exceeds the numerical rank of the supplied data sketch.');
end
basis=modes(:,1:target_rank);

% Measure the projection residual block by block without a second data copy
residual_sq=0;
for n=1:numel(blocks)
    residual_sq=residual_sq+norm(blocks{n}-basis*(basis'*blocks{n}),'fro')^2;
end
result=struct('basis',basis,'singular_vals',singular_vals,...
              'total_energy',total_energy,'residual_abs',sqrt(residual_sq),...
              'residual_rel',sqrt(residual_sq/total_energy),...
              'orth_error',norm(basis'*basis-eye(target_rank),'fro'));

end

% Consistency enforcement
function grumble(blocks,weights,target_rank,options)

if ~iscell(blocks)||~iscolumn(blocks)||isempty(blocks)
    error('blocks must be a nonempty cell column.');
end
for n=1:numel(blocks)
    if ~isa(blocks{n},'double')||~ismatrix(blocks{n})||isempty(blocks{n})||...
       ~all(isfinite(blocks{n}),'all')||size(blocks{n},1)~=size(blocks{1},1)
        error('blocks must contain finite nonempty double matrices with common row count.');
    end
end
if ~isa(weights,'double')||~isreal(weights)||~isequal(size(weights),size(blocks))||...
   ~all(isfinite(weights))||any(weights<0)||~any(weights>0)
    error('weights must be a matching finite nonnegative double column with a positive entry.');
end
ndims_data=min(size(blocks{1},1),sum(cellfun(@(x)size(x,2),blocks)));
if ~isnumeric(target_rank)||~isreal(target_rank)||~isscalar(target_rank)||...
   ~isfinite(target_rank)||target_rank<1||mod(target_rank,1)~=0||target_rank>ndims_data
    error('target_rank must be a positive integer within the data dimensions.');
end
fields={'method','oversampling','power_iters','seed'};
if ~isstruct(options)||~isscalar(options)||~all(isfield(options,fields))
    error('options must supply method, oversampling, power_iters, and seed.');
end
if ~ischar(options.method)||~isrow(options.method)||...
   ~ismember(options.method,{'svd','randomised'})
    error('method must be svd or randomised.');
end
for n=2:4
    value=options.(fields{n});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||...
       value<0||mod(value,1)~=0
        error('oversampling, power_iters, and seed must be nonnegative integers.');
    end
end
if options.seed>2^32-1
    error('seed must not exceed 2^32-1.');
end
if strcmp(options.method,'svd')&&(options.oversampling~=0||options.power_iters~=0)
    error('svd requires zero oversampling and power_iters.');
end
if strcmp(options.method,'randomised')&&target_rank+options.oversampling>ndims_data
    error('target_rank plus oversampling exceeds the data dimensions.');
end

end


