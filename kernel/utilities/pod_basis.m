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
% caller's random state is preserved, including on an error. Randomised
% construction applies weighted blocks without a full snapshot concatenation;
% workspace includes n-by-sketch_rank and active_columns-by-sketch_rank arrays
% for range products, compressed data, and residual measurement. Zero-weight
% blocks are validated, but excluded from construction. Exact singular-value
% vectors retain min(n,total supplied columns) entries; randomised vectors
% retain target_rank+oversampling entries, padding excluded directions with
% zeros. Internal sketches are limited by the positively weighted columns.
% Numerical rank tolerance still uses the original supplied-data dimensions.
% Randomised sketches remain real when all positively weighted blocks are real.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pod_basis.m>

function result=pod_basis(blocks,weights,target_rank,options)

% Check the supplied data contract
grumble(blocks,weights,target_rank,options);

% Exclude zero-weight blocks while retaining the supplied-data dimensions
row_count=size(blocks{1},1); data_columns=sum(cellfun(@(x)size(x,2),blocks));
keep_blocks=weights>0; blocks=blocks(keep_blocks); weights=weights(keep_blocks);
column_count=sum(cellfun(@(x)size(x,2),blocks));

% Measure the weighted data norm without concatenating or copying the blocks
data_norm=0; block_width=min(target_rank+options.oversampling,column_count);
for n=1:numel(blocks)
    for k=1:block_width:size(blocks{n},2)
        columns=k:min(k+block_width-1,size(blocks{n},2));
        data_norm=hypot(data_norm,norm(sqrt(weights(n))*blocks{n}(:,columns),'fro'));
    end
end
total_energy=data_norm^2;
if ~isfinite(total_energy)||total_energy==0
    error('weighted data must have finite nonzero energy.');
end

% Compute the exact basis or a seeded blockwise randomised range sketch
if strcmp(options.method,'svd')
    output_size=min(row_count,data_columns);
    data=zeros(row_count,column_count); first_column=1;
    for n=1:numel(blocks)
        columns=first_column:first_column+size(blocks{n},2)-1;
        data(:,columns)=sqrt(weights(n))*blocks{n};
        first_column=first_column+size(blocks{n},2);
    end
    [modes,values]=svd(data,'econ');
else
    output_size=target_rank+options.oversampling;
    saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng));
    rng(options.seed,'twister'); sketch_rank=block_width;
    omega=randn(column_count,sketch_rank);
    if ~all(cellfun(@isreal,blocks))
        omega=complex(omega,randn(column_count,sketch_rank));
    end
    for k=0:options.power_iters

        % Apply the weighted snapshot matrix using only the supplied blocks
        range=zeros(row_count,sketch_rank); first_column=1;
        for n=1:numel(blocks)
            columns=first_column:first_column+size(blocks{n},2)-1;
            range=range+blocks{n}*(sqrt(weights(n))*omega(columns,:));
            first_column=first_column+size(blocks{n},2);
        end
        [range,~]=qr(range,0);
        if k<options.power_iters

            % Reorthogonalise the adjoint range before the next power iteration
            omega=zeros(column_count,sketch_rank); first_column=1;
            for n=1:numel(blocks)
                columns=first_column:first_column+size(blocks{n},2)-1;
                omega(columns,:)=blocks{n}'*(sqrt(weights(n))*range);
                first_column=first_column+size(blocks{n},2);
            end
            [omega,~]=qr(omega,0);

        end

    end

    % Factor only the compressed snapshots, without right singular vectors
    data=zeros(sketch_rank,column_count); first_column=1;
    for n=1:numel(blocks)
        columns=first_column:first_column+size(blocks{n},2)-1;
        data(:,columns)=(sqrt(weights(n))*range')*blocks{n};
        first_column=first_column+size(blocks{n},2);
    end
    [left_modes,values]=svd(data,'econ'); modes=range*left_modes;
    clear restore_rng
end

% Reject unsupported rank rather than returning arbitrary null-space modes
singular_vals=[diag(values);zeros(output_size-size(values,1),1)];
rank_tol=max(row_count,data_columns)*eps(singular_vals(1));
if singular_vals(target_rank)<=rank_tol
    error('target_rank exceeds the numerical rank of the supplied data sketch.');
end
basis=modes(:,1:target_rank);

% Measure the projection residual in column slices bounded by the sketch rank
residual_abs=0;
for n=1:numel(blocks)
    for k=1:block_width:size(blocks{n},2)
        columns=k:min(k+block_width-1,size(blocks{n},2));
        data=sqrt(weights(n))*blocks{n}(:,columns);
        residual_abs=hypot(residual_abs,norm(data-basis*(basis'*data),'fro'));
    end
end
result=struct('basis',basis,'singular_vals',singular_vals,...
              'total_energy',total_energy,'residual_abs',residual_abs,...
              'residual_rel',residual_abs/data_norm,...
              'orth_error',norm(basis'*basis-eye(target_rank),'fro'));

end

% Consistency enforcement
function grumble(blocks,weights,target_rank,options)
if ~iscell(blocks)||~iscolumn(blocks)||isempty(blocks)
    error('blocks must be a nonempty cell column.');
end
for n=1:numel(blocks)
    if ~isa(blocks{n},'double')||~ismatrix(blocks{n})||isempty(blocks{n})||...
       any(isnan(blocks{n}),'all')||any(isinf(blocks{n}),'all')||...
       size(blocks{n},1)~=size(blocks{1},1)
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


