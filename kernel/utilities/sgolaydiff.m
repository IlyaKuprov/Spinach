% Savitzky-Golay differentiation of noisy sampled signals by local
% least-squares polynomial fitting. Syntax:
%
%          dy=sgolaydiff(y,der_order,npoints,poly_order)
%
% Parameters:
%
%    y          - N-by-M signal matrix; rows are samples and
%                 columns are independent signals
%
%    der_order  - derivative order; order 0 returns the
%                 smoothed signal
%
%    npoints    - odd number of points in the local least-
%                 squares window
%
%    poly_order - order of the local polynomial
%
% Outputs:
%
%    dy         - N-by-M derivative matrix on a unit-step
%                 uniform grid
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sgolaydiff.m>

function dy=sgolaydiff(y,der_order,npoints,poly_order)

% Check consistency
grumble(y,der_order,npoints,poly_order);

% Count the samples
nsamps=size(y,1);

% Preallocate the derivative estimate
dy=0*y;

% Compute sided and centred local fits
win_half=(npoints-1)/2; fac=factorial(der_order);
for n=1:nsamps

    % Select a local window around the current sample
    win_left=max(1,min(n-win_half,nsamps-npoints+1));
    win_idx=win_left:(win_left+npoints-1);

    % Centre and scale the local sample offsets
    loc_grid=win_idx(:)-n;
    grid_scale=max(abs(loc_grid));
    loc_grid=loc_grid/grid_scale;

    % Build the local Vandermonde matrix
    V=ones(npoints,poly_order+1);
    for k=1:poly_order
        V(:,k+1)=loc_grid.^k;
    end

    % Fit the local polynomial by QR-backed least squares
    coeff=V\y(win_idx,:);

    % Extract the requested derivative at the expansion point
    dy(n,:)=fac*coeff(der_order+1,:)/(grid_scale^der_order);

end

end

% Consistency enforcement
function grumble(y,der_order,npoints,poly_order)
if (~isfloat(y))||isempty(y)||issparse(y)||(~ismatrix(y))
    error('y must be a non-empty dense floating-point matrix.');
end
if any(~isfinite(y(:)))
    error('y must not contain NaN or Inf values.');
end
if (~isnumeric(der_order))||(~isreal(der_order))||(numel(der_order)~=1)||...
   (der_order<0)||(mod(der_order,1)~=0)
    error('der_order must be a non-negative real integer.');
end
if (~isnumeric(npoints))||(~isreal(npoints))||(numel(npoints)~=1)||...
   (npoints<3)||(mod(npoints,1)~=0)
    error('npoints must be a positive odd integer greater than two.');
end
if mod(npoints,2)~=1
    error('npoints must be an odd integer.');
end
if (~isnumeric(poly_order))||(~isreal(poly_order))||(numel(poly_order)~=1)||...
   (poly_order<0)||(mod(poly_order,1)~=0)
    error('poly_order must be a non-negative real integer.');
end
nsamps=size(y,1);
if nsamps<3
    error('y must contain at least three sample rows.');
end
if npoints>nsamps
    error('npoints must not exceed the number of samples.');
end
if der_order>poly_order
    error('der_order must not exceed poly_order.');
end
if poly_order>=npoints
    error('poly_order must be smaller than npoints.');
end
end

% There are only two hard things in Computer Science: cache
% invalidation and naming things.
%
% Phil Karlton

