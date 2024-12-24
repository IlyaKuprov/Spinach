% Tikhonov regularised solution to K*x=y with a positivity const-
% raint on x using regularised Newton-Raphson method. Syntax:
%
%          [x,err,reg]=tikhonov(K,D,KtK,DtD,H,y,lambda)
%
% Parameters:
%
%    K      - kernel matrix, may be complex, may be non-square
%
%    D      - regularisation matrix, leave empty to use finite
%             difference second derivative matrix
%
%    KtK    - K'*K, for repeated calls it may be faster to pre-
%             compute this quantity, leave empty otherwise
%
%    DtD    - D'*D, for repeated calls it may be faster to pre-
%             compute this quantity, leave empty otherwise
%
%    H      - Tikhonov Hessian 2*real(KtK+lambda*DtD), for re-
%             peated calls it may be faster to precompute this
%             quantity, leave empty otherwise
%
%    y      - a column vector, may be complex
%
%    lambda - Tikhonov regularisation parameter
%
% Outputs:
%
%    x      - a real vector, a minimum (subject to positivity)
%             of norm(K*x-y,2)^2+lambda*norm(D*x,2)^2
%
%    err    - error signal norm(K*x-y,2)^2
%
%    reg    - regularisation signal norm(D*x,2)^2
%
% Note: for best numerical performance, scale K to have approxima-
%       tely unit 2-norm, and y to have approximately unit 1-norm.
%
% Note: see tikhoind.m for the indeterminate solver.
%
% ilya.kuprov@weizmann.ac.uk
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=tikhonov.m>

function [x,err,reg]=tikhonov(K,D,KtK,DtD,H,y,lambda)

% Check consistency
grumble(K,D,KtK,DtD,H,y,lambda);

% Find output size
outsize=size(K,2);

% Default D is the second derivative matrix
if isempty(D), D=fdmat(outsize,5,2,'wall'); end

% Various repeating objects
if isempty(KtK), KtK=K'*K; end
if isempty(DtD), DtD=D'*D; end
if isempty(H), H=2*real(KtK+lambda*DtD); end

    % Tikhonov regularised error signal
    function [tikh,grad]=lsq_err(x,y)
        
        % Composite Tikhonov error signal
        tikh=norm(K*x-y,2)^2+lambda*norm(D*x,2)^2;

        % Gradient only computed if requested
        if nargout>1, grad=2*real(KtK*x-K'*y+lambda*(DtD*x)); end
       
    end

% Optimisation options and bounds
options=optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
                     'MaxIterations',100,'MaxFunctionEvaluations',inf,...
                     'SpecifyObjectiveGradient',true,'HessianFcn',@(varargin)H);
lbound=zeros(outsize,1); ubound=inf(outsize,1);

% Run the optimisation
x=fmincon(@(x)lsq_err(x,y),ones(outsize,1),[],[],[],[],lbound,ubound,[],options);

% Error and regularisation signals
if nargout>1, err=norm(K*x-y,2)^2; end
if nargout>2, reg=norm(D*x,2)^2; end

end

% Consistency enforcement
function grumble(K,D,KtK,DtD,H,y,lambda)
if (~isnumeric(K))||(~isnumeric(D))||(~isnumeric(KtK))||...
   (~isnumeric(DtD))||(~isnumeric(H))||(~isnumeric(y))||...
   (~isnumeric(lambda))
    error('all inputs must be numeric.');
end
if size(K,1)~=size(y,1)
    error('dimensions of K and y are not consistent.');
end
if ~isempty(KtK)
    if size(KtK,1)~=size(KtK,2)
        error('KtK must be a square matrix.');
    end
end
if ~isempty(DtD)
    if size(DtD,1)~=size(DtD,2)
        error('DtD must be a square matrix.');
    end
end
if ~isempty(H)
    if size(H,1)~=size(H,2)
        error('H must be a square matrix.');
    end
end
if (~isreal(lambda))||(~isscalar(lambda))||(lambda<0)
    error('lambda must be a positive real scalar.');
end
end

% All generalisations are untrue, including this one.

