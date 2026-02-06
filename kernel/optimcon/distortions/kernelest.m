% FIR convolution kernel estimation from input/output samples.
%
%          h=kernelest(x,y,ker_len,method,align,lambda)
%
% Parameters:
%
%    x        - input samples on a uniform grid
%
%    y        - output samples on the same grid
%
%    ker_len  - kernel length (number of taps)
%
%    method   - 'backslash' (default) | 'pinv' | 'svd' | 'tikh'
%
%    align    - 'causal' (default) or 'same' output alignment
%
%    lambda   - Tikhonov parameter for 'tikh' (optional)
%
% Outputs:
%
%    h        - estimated convolution kernel
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kernelest.m>

function h=kernelest(x,y,ker_len,method,align,lambda)

% Set the defaults
if nargin<4, method='backslash'; end
if nargin<5, align='causal'; end
if nargin<6, lambda=1e-6; end

% Check consistency
grumble(x,y,ker_len,method,align,lambda);

% Column vectors and sample count
x=x(:); y=y(:); npts=numel(x);

% Build convolution matrix
conv_mat=toeplitz([x;   zeros(ker_len-1,1)],...
                  [x(1) zeros(1,ker_len-1)]);

% Alignment choice
switch lower(align)

    case 'causal'

        % Select the first npts rows
        sys_mat=conv_mat(1:npts,:);

    case 'same'

        % Select the central npts rows
        row_start=floor(ker_len/2)+1;
        sys_mat=conv_mat(row_start:(row_start+npts-1),:);

    otherwise

        % Complain and bomb out
        error('align must be ''causal'' or ''same''.');

end

% Compute kernel 
switch lower(method)

    case 'backslash'

        % LLSq
        h=sys_mat\y;

    case 'pinv'

        % Pseudoinverse
        h=pinv(sys_mat)*y;

    case 'svd'

        % Truncated-SVD pseudoinverse
        [u_mat,s_mat,v_mat]=svd(sys_mat,'econ');
        s_vals=diag(s_mat);
        s_tol=max(size(sys_mat))*eps(max(s_vals));
        s_inv=zeros(size(s_vals));
        s_inv(s_vals>s_tol)=1./s_vals(s_vals>s_tol);
        h=v_mat*(s_inv.*(u_mat'*y));

    case 'tikh'

        % Solve Tikhonov-regularised system
        h=(sys_mat'*sys_mat+lambda*eye(ker_len))\(sys_mat'*y);

    otherwise

        % Complain and bomb out
        error('Unknown method: %s',method);

end

end

% Consistency enforcement
function grumble(x,y,ker_len,method,align,lambda)
if (~isnumeric(x))||(~isvector(x))
    error('x must be a numeric vector.');
end
if (~isnumeric(y))||(~isvector(y))||(numel(y)~=numel(x))
    error('y must be a numeric vector with the same length as x.');
end
if (~isnumeric(ker_len))||(~isreal(ker_len))||...
   (~isscalar(ker_len))||(mod(ker_len,1)~=0)||(ker_len<1)
    error('ker_len must be a positive integer.');
end
if (~ischar(method))&&(~isstring(method))
    error('method must be a character string.');
end
if (~ischar(align))&&(~isstring(align))
    error('align must be a character string.');
end
if (~isnumeric(lambda))||(~isreal(lambda))||...
   (~isscalar(lambda))||(lambda<=0)
    error('lambda must be a positive real scalar.');
end
end

% The recipe is formalised here, for the sake of
% cookery-column propriety, but please know that
% this is a pie which lends itself to flexibili-
% to and approximation.
%
% Olivia Potts, the cooking
% columninst at The Spectator

