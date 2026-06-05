% Reconstructs the instantaneous Bloch equation rotation axis of
% from a 3D magnetisation trajectory. Syntax:
%
%                 [ax,ay,az]=bloch_axis(x,y,z)
%
% Parameters:
%
%   x, y, z    - row vectors of equal length con-
%                taining the trajectory
%
% Output:
%
%   ax, ay, az - row vectors of equal length con-
%                taining the instantaneous axis 
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bloch_axis.m>

function [ax,ay,az]=bloch_axis(x,y,z)

% Check consistency
grumble(x,y,z);

% Get first and second derivatives
dx_dt=fdvec(x,5,1); d2x_dt2=fdvec(x,5,2);
dy_dt=fdvec(y,5,1); d2y_dt2=fdvec(y,5,2);
dz_dt=fdvec(z,5,1); d2z_dt2=fdvec(z,5,2);

% Get instantaneous rotation axis
ax=dz_dt.*d2y_dt2-dy_dt.*d2z_dt2;
ay=dx_dt.*d2z_dt2-dz_dt.*d2x_dt2;
az=dy_dt.*d2x_dt2-dx_dt.*d2y_dt2;

end

% Consistency enforcement
function grumble(x,y,z)
if (~isnumeric(x))||(~isreal(x))||...
   (~ismatrix(x))||any(~isfinite(x),'all')
    error('x must be a finite real matrix.');
end
if (~isnumeric(y))||(~isreal(y))||...
   (~ismatrix(y))||any(~isfinite(y),'all')
    error('y must be a finite real matrix.');
end
if (~isnumeric(z))||(~isreal(z))||...
   (~ismatrix(z))||any(~isfinite(z),'all')
    error('z must be a finite real matrix.');
end
if (~isequal(size(x),size(y)))||...
   (~isequal(size(x),size(z)))
    error('x, y, and z must have the same dimensions.');
end
end

% The Proton VPN servers for Afghanistan running at 
% over 80% load was not something that we expected
% to see when we set them up last year. 
%
% David Peterson, after the UK started
% requiring proof of age on porn sites
% in August 2025

