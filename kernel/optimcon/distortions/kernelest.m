%RECOVER_CONV_KERNEL  Estimate an L-tap FIR convolution kernel h from x and y.
%
% Model (default, causal):   y(n) ≈ sum_{k=1..L} h(k) * x(n-k+1)   (zero-padded for n-k+1<=0)
%
% Inputs
%   x,y    : input/output samples on the same grid (same length N)
%   L      : kernel length (number of taps / grid points)
%   method : 'backslash' (default) | 'pinv' | 'svd' | 'tikh'
%   align  : 'causal' (default) or 'same' (matches conv(x,h,'same') indexing)
%   lambda : Tikhonov regularization for 'tikh' (optional, e.g. 1e-6)
%
% Output
%   h      : estimated kernel (Lx1)

function h = kernelest(x,y,L,method,align,lambda)

x = x(:);  y = y(:);
N = numel(x);

if nargin < 4 || isempty(method), method = 'backslash'; end
if nargin < 5 || isempty(align),  align  = 'causal';    end

% Full convolution matrix: conv(x,h) = Xfull*h, size (N+L-1) x L
Xfull = toeplitz([x; zeros(L-1,1)], [x(1) zeros(1,L-1)]);

% Pick the rows that correspond to how your instrument outputs "same-grid" data
switch lower(align)
    case 'causal'
        A = Xfull(1:N, :);                         % y ≈ first N samples of conv(x,h)
    case 'same'
        start = floor(L/2) + 1;                    % typical conv(...,'same') indexing
        A = Xfull(start:start+N-1, :);             % y ≈ central N samples of conv(x,h)
    otherwise
        error('align must be ''causal'' or ''same''.');
end

% Solve least-squares for h
switch lower(method)
    case 'backslash'
        h = A \ y;                                 % LS fit (or exact if square)
    case 'pinv'
        h = pinv(A) * y;                           % Moore-Penrose pseudoinverse
    case 'svd'
        [U,S,V] = svd(A,'econ');                   % truncated-SVD pseudoinverse
        s = diag(S);
        tol = max(size(A)) * eps(max(s));
        sInv = zeros(size(s));
        sInv(s > tol) = 1 ./ s(s > tol);
        h = V * (sInv .* (U' * y));
    case 'tikh'
        if nargin < 6 || isempty(lambda), lambda = 1e-6; end
        h = (A'*A + lambda*eye(L)) \ (A'*y);        % ridge / Tikhonov
    otherwise
        error('Unknown method: %s', method);
end
end

% The recipe is formalised here, for the sake of
% cookery-column propriety, but please know that
% this is a pie which lends itself to flexibili-
% to and approximation.
%
% Olivia Potts, the cooking
% columninst at The Spectator

