% Non-orthogonal channel distortion model. Treats odd rows of
% multi-row waveform arrays as in-phase channels, and even rows
% as out-of-phase channels. The in-phase channel is kept fixed;
% the out-of-phase channel is tilted so that its true angle to
% the in-phase channel is user-specified. Syntax:
%
%                    [w,J]=non_orth(w,xy_ang)
%
% Parameters:
%
%    w        - waveform, one time slice per column, and
%               rows arranged as XYXY... with respect to
%               in-phase and quadrature parts on each
%               control channel
%
%    xy_ang   - true angle, in degrees, between the instru-
%               ment output directions of each X,Y control
%               pair; may be a scalar or one value per pair,
%               with 90 degrees corresponding to no distortion
%
% Outputs:
%
%    w        - distorted waveform, same dimension as the
%               input waveform
%
%    J        - Jacobian matrix with respect to vectorisa-
%               tions of the output and the input arrays
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=non_orth.m>

function [w,J]=non_orth(w,xy_ang)

% Check consistency
grumble(w,xy_ang);

% Count waveform dimensions
[nrows,ncols]=size(w); nchannels=nrows/2;

% Expand a scalar angle specification
if isscalar(xy_ang), xy_ang=xy_ang*ones(1,nchannels); end

% Preallocate the output waveform
w_dist=zeros(size(w),'like',w);

% Preallocate the channel mixing matrix triplets
if nargout>1
    row_idx=zeros(1,3*nchannels);
    col_idx=zeros(1,3*nchannels);
    mat_val=zeros(1,3*nchannels);
end

% Loop over control channel pairs
for n=1:nchannels

    % Get the row numbers
    x_row=2*n-1; y_row=2*n;

    % Get the trigonometric factors
    cos_ang=cosd(xy_ang(n)); sin_ang=sind(xy_ang(n));

    % Apply the non-orthogonal channel tilt
    w_dist(x_row,:)=w(x_row,:)+cos_ang*w(y_row,:);
    w_dist(y_row,:)=sin_ang*w(y_row,:);

    % Store the Jacobian block triplets
    if nargout>1
        elem_idx=3*n-2;
        row_idx(elem_idx:(elem_idx+2))=[x_row x_row y_row];
        col_idx(elem_idx:(elem_idx+2))=[x_row y_row y_row];
        mat_val(elem_idx:(elem_idx+2))=[1 cos_ang sin_ang];
    end

end

% Return distorted waveform
w=w_dist;

% Build the vectorised Jacobian
if nargout>1
    chan_mat=sparse(row_idx,col_idx,mat_val,nrows,nrows);
    J=kron(speye(ncols),chan_mat);
end

end

% Consistency enforcement
function grumble(w,xy_ang)
if (~isnumeric(w))||(~isreal(w))
    error('w must be an array of real numbers.');
end
if mod(size(w,1),2)~=0
    error('the number of rows in w must be even.');
end
if (~isnumeric(xy_ang))||(~isreal(xy_ang))||isempty(xy_ang)||...
   ((numel(xy_ang)~=1)&&(numel(xy_ang)~=size(w,1)/2))
    error('xy_ang must be a real scalar, or have one element per X,Y control pair.');
end
if any(~isfinite(xy_ang),'all')||any(xy_ang<=0,'all')||any(xy_ang>=180,'all')
    error('xy_ang must contain angles strictly between 0 and 180 degrees.');
end
end


