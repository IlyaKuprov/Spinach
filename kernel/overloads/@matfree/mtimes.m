% Scale a matrix-free core or apply it to a full block of columns.
% Syntax:
%
%                         answer=mtimes(left,right)
%
% Parameters:
%
%    left,right - a matfree core and a finite numeric scalar, or a
%                 matfree core followed by a full columns-by-n block
%
% Outputs:
%
%    answer     - scaled core or full rows-by-n block
%
% ilya.kuprov@weizmann.ac.il

function answer=mtimes(left,right)

% Check consistency
grumble(left,right);

% Scale an implicit operator from either side
if isnumeric(left)&&~isa(left,'matfree')&&isscalar(left)
    answer=right; answer.coeff=left*right.coeff;
elseif isnumeric(right)&&~isa(right,'matfree')&&isscalar(right)
    answer=left; answer.coeff=left.coeff*right;
else

    % Apply the forward action to every right-hand side
    answer=left.coeff*left.forward(right);
end

end

% Consistency enforcement
function grumble(left,right)
if isa(left,'matfree')&&isnumeric(right)&&~isa(right,'matfree')
    if isscalar(right)
        if ~isfinite(right), error('scalar coefficient must be finite.'); end
    elseif ~ismatrix(right)||issparse(right)||(size(right,1)~=left.dims(2))
        error('right block must be full with one row per operator column.');
    end
elseif isa(right,'matfree')&&isnumeric(left)&&isscalar(left)
    if ~isfinite(left), error('scalar coefficient must be finite.'); end
else
    error('use a numeric scalar or a full block with a matrix-free core.');
end
end


