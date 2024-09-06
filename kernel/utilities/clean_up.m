% Array clean-up utility. Drops non-zero elements with magnitude below the
% user-specified tolerance and converts between sparse and full storage de-
% pending on the density of non-zeroes in the array. Syntax:
%
%                 A=clean_up(spin_system,A,nonzero_tol)
%
% Parameters:
%
%            A  - a numerical array or a cell array thereof
%
%   nonzero_tol - nonzero tolerance
%
% Outputs:
%
%            A  - cleaned-up array
%
% i.kuprov@soton.ac.uk
% alex_nevzorov@ncsu.edu
%
% <https://spindynamics.org/wiki/index.php?title=clean_up.m>

function A=clean_up(spin_system,A,nonzero_tol)

% Skip opium objects
if isa(A,'opium'), return; end

% Skip if disabled
if (nonzero_tol==0)||isnan(nonzero_tol), return; end

% Process cells recursively
if iscell(A)
    for n=1:numel(A)
        A{n}=clean_up(spin_system,A{n},nonzero_tol);
    end
    return
end

% Process polyadics recursively
if isa(A,'polyadic')
    for n=1:numel(A.prefix)
        A.prefix{n}=clean_up(spin_system,A.prefix{n},nonzero_tol);
    end
    for n=1:numel(A.suffix)
        A.suffix{n}=clean_up(spin_system,A.suffix{n},nonzero_tol);
    end
    for n=1:numel(A.cores)
        for k=1:numel(A.cores{n})
            A.cores{n}{k}=clean_up(spin_system,A.cores{n}{k},nonzero_tol);
        end
    end
    return
end
   
% Check consistency
grumble(A,nonzero_tol);

% Check if clean-up is allowed
if ~ismember('clean-up',spin_system.sys.disable)

    % Use MEX if appropriate
    if ismember('mex',spin_system.sys.enable)&&...
       issparse(A)&&(~isa(A,'gpuArray'))

        % Memory-friendly in-place MEX
        prune_cpu(A,nonzero_tol); A=1*A;

    else

        % Use the memory-hungry generic method
        A=nonzero_tol*round((1/nonzero_tol)*A);

    end
    
    % A small non-zero matrix should always be full
    if issparse(A)&&(nnz(A)>0)&&...
       any(size(A)<spin_system.tols.small_matrix), A=full(A); end
    
    % A big sparse matrix with too many non-zeros should be full
    if issparse(A)&&(nnz(A)/numel(A)>spin_system.tols.dense_matrix), A=full(A); end
    
    % A big full matrix with too few non-zeros should be sparse
    if (~issparse(A))&&(nnz(A)/numel(A)<spin_system.tols.dense_matrix)&&...
                       (all(size(A)>spin_system.tols.small_matrix))
        A=sparse(A);
    end
    
end

end

% Consistency enforcement
function grumble(A,nonzero_tol)
if (~isnumeric(A))
    error('A must be numeric.');
end
if ~isnumeric(nonzero_tol)
    error('the tolerance parameter must be numeric.');
end
if (~isreal(nonzero_tol))||(numel(nonzero_tol)~=1)||(nonzero_tol<=0)
    error('nonzero_tol parameter must be a positive real number.');
end
end

% The most dangerous man to any government is the man who is able to think
% things out for himself, without regard to the prevailing superstitions
% and taboos.
%
% H.L. Mencken

