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
%     result.total_energy  - squared Frobenius norm of weighted supplied data;
%                            may round to zero or infinity
%     result.residual_abs  - ||data-basis*(basis'*data)||_F, on supplied data
%     result.residual_rel  - residual_abs divided by the unsquared data norm
%     result.orth_error    - ||basis'*basis-eye(target_rank)||_F
%     result.scale_exponent - common binary exponent of weighted data
%     result.singular_vals_scaled - singular values before binary rescaling
%     result.data_norm_scaled - Frobenius norm before binary rescaling
%
% Exact SVD minimises the supplied weighted snapshot Frobenius projection
% error. Randomised sketches and pre-truncated factors qualify that claim.
% Residuals are measured data diagnostics, not observable error bounds.
% Rank-deficient inputs are rejected if the requested rank exceeds the
% numerical rank; no arbitrary null-space directions are supplied. The
% caller's random state is preserved, including on an error. Factorisation,
% rank and relative residuals use a power-of-two-normalised weighted matrix.
% Squared physical energy never decides rank. Physical diagnostics may exceed
% double range; the scaled diagnostics remain available. Loss of a nonzero
% component during normalisation is rejected explicitly. Randomised
% construction applies weighted blocks without a full snapshot concatenation;
% workspace includes n-by-sketch_rank and active_columns-by-sketch_rank arrays
% for range products, compressed data, and residual measurement. Zero-weight
% blocks are validated, but excluded from construction. Exact singular-value
% vectors retain min(n,total supplied columns) entries; randomised vectors
% retain target_rank+oversampling entries, padding excluded directions with
% zeros. Internal sketches are limited by the positively weighted columns.
% Numerical rank tolerance still uses the original supplied-data dimensions.
% Sketch coefficients are released before compressed snapshots are allocated.
% Only retained randomised modes are formed, after numerical rank validation.
% Randomised sketches remain real when all positively weighted blocks are real.
%
% Algorithms: Benner et al., DOI:10.1137/130932715; Halko et al.,
% DOI:10.1137/090771806, Algorithms 4.4 and 5.1. All integer inputs are double.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pod_basis.m>

function result=pod_basis(blocks,weights,target_rank,options)

% Check the supplied data contract
grumble(blocks,weights,target_rank,options);

% Retain original dimensions but exclude zero-weight blocks from the operator
row_count=size(blocks{1},1); data_columns=sum(cellfun(@(x)size(x,2),blocks));
keep=weights>0; blocks=blocks(keep); weights=weights(keep);
column_count=sum(cellfun(@(x)size(x,2),blocks));
sketch_rank=min(target_rank+options.oversampling,column_count);

% Represent weighted block scales without forming their potentially large product
[weight_frac,weight_exp]=log2(sqrt(weights)); block_exp=zeros(size(weights));
for n=1:numel(blocks)
    largest=0;
    for k=1:sketch_rank:size(blocks{n},2)
        columns=k:min(k+sketch_rank-1,size(blocks{n},2)); data=blocks{n}(:,columns);
        largest=max(largest,max(max(abs(real(data)),[],'all'),max(abs(imag(data)),[],'all')));
    end
    if largest==0, block_exp(n)=-Inf; else, [~,block_exp(n)]=log2(largest); end
end
clear data
scale_exponent=max(block_exp+weight_exp);
if ~isfinite(scale_exponent), error('weighted data must be nonzero.'); end

% Measure the normalised operator norm in bounded column slices
data_norm=0;
for n=1:numel(blocks)
    for k=1:sketch_rank:size(blocks{n},2)
        columns=k:min(k+sketch_rank-1,size(blocks{n},2));
        data=weighted_slice(blocks{n}(:,columns),weight_frac(n),weight_exp(n)-scale_exponent);
        data_norm=hypot(data_norm,norm(data,'fro'));
    end
end

% Factor the complete normalised data or construct its stabilised range sketch
if strcmp(options.method,'svd')
    output_size=min(row_count,data_columns); data=zeros(row_count,column_count);
    first_column=1;
    for n=1:numel(blocks)
        for k=1:sketch_rank:size(blocks{n},2)
            columns=k:min(k+sketch_rank-1,size(blocks{n},2));
            data(:,first_column+columns-1)=weighted_slice(blocks{n}(:,columns),...
                                                         weight_frac(n),weight_exp(n)-scale_exponent);
        end
        first_column=first_column+size(blocks{n},2);
    end
    [modes,singular_vals]=svd(data,'econ','vector'); clear data
else
    output_size=target_rank+options.oversampling;
    saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng)); rng(options.seed,'twister');
    omega=randn(column_count,sketch_rank);
    if ~all(cellfun(@isreal,blocks)), omega=complex(omega,randn(column_count,sketch_rank)); end
    for j=0:options.power_iters

        % Apply the normalised weighted data and reorthogonalise its range
        Q=zeros(row_count,sketch_rank); first_column=1;
        for n=1:numel(blocks)
            for k=1:sketch_rank:size(blocks{n},2)
                columns=k:min(k+sketch_rank-1,size(blocks{n},2));
                data=weighted_slice(blocks{n}(:,columns),weight_frac(n),weight_exp(n)-scale_exponent);
                Q=Q+data*omega(first_column+columns-1,:);
            end
            first_column=first_column+size(blocks{n},2);
        end
        [Q,~]=qr(Q,0); clear data
        if j<options.power_iters

            % Apply the adjoint and reorthogonalise before another data pass
            first_column=1;
            for n=1:numel(blocks)
                for k=1:sketch_rank:size(blocks{n},2)
                    columns=k:min(k+sketch_rank-1,size(blocks{n},2));
                    data=weighted_slice(blocks{n}(:,columns),weight_frac(n),weight_exp(n)-scale_exponent);
                    omega(first_column+columns-1,:)=data'*Q;
                end
                first_column=first_column+size(blocks{n},2);
            end
            [omega,~]=qr(omega,0); clear data

        end

    end

    % Compress the normalised data after releasing Gaussian coefficients
    clear omega
    data=zeros(sketch_rank,column_count); first_column=1;
    for n=1:numel(blocks)
        for k=1:sketch_rank:size(blocks{n},2)
            columns=k:min(k+sketch_rank-1,size(blocks{n},2));
            data(:,first_column+columns-1)=Q'*weighted_slice(blocks{n}(:,columns),...
                                                           weight_frac(n),weight_exp(n)-scale_exponent);
        end
        first_column=first_column+size(blocks{n},2);
    end
    [modes,singular_vals]=svd(data,'econ','vector'); clear data restore_rng
end

% Apply the original-dimension rank threshold in normalised coordinates
singular_vals=[singular_vals;zeros(output_size-numel(singular_vals),1)];
rank_tol=max(row_count,data_columns)*eps(singular_vals(1));
if target_rank>numel(singular_vals)||singular_vals(target_rank)<=rank_tol
    error('target_rank exceeds the numerical rank of the supplied data sketch.');
end
basis=modes(:,1:target_rank);
if strcmp(options.method,'randomised'), basis=Q*basis; end
clear modes Q

% Measure the actual supplied-data residual rather than a spectral tail estimate
residual=0;
for n=1:numel(blocks)
    for k=1:sketch_rank:size(blocks{n},2)
        columns=k:min(k+sketch_rank-1,size(blocks{n},2));
        data=weighted_slice(blocks{n}(:,columns),weight_frac(n),weight_exp(n)-scale_exponent);
        residual=hypot(residual,norm(data-basis*(basis'*data),'fro'));
    end
end
result=struct('basis',basis,'singular_vals',pow2(singular_vals,scale_exponent),...
              'total_energy',pow2(data_norm^2,2*scale_exponent),...
              'residual_abs',pow2(residual,scale_exponent),'residual_rel',residual/data_norm,...
              'orth_error',norm(basis'*basis-speye(target_rank),'fro'),...
              'scale_exponent',scale_exponent,'singular_vals_scaled',singular_vals,...
              'data_norm_scaled',data_norm);

end

% Apply one common binary scale to a weighted data slice without losing components
function data=weighted_slice(data,fraction,exponent)
real_data=pow2(real(data),exponent)*fraction;
imag_data=pow2(imag(data),exponent)*fraction;
if any((real(data)~=0)&(real_data==0),'all')||any((imag(data)~=0)&(imag_data==0),'all')
    error('weighted data dynamic range exceeds the normalised double representation.');
end
if isreal(data), data=real_data; else, data=complex(real_data,imag_data); end
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
if ~isa(target_rank,'double')||~isreal(target_rank)||~isscalar(target_rank)||...
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
    if ~isa(value,'double')||~isreal(value)||~isscalar(value)||~isfinite(value)||...
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


