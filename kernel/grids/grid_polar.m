% Generates a balanced polar grid in which the density of 
% points does not increase towards the centre. Syntax:
%
%          [phi,r,L]=grid_polar(ncircles,rmax)
%
% Parameters:
%
%     ncircles - number of radial circles 
%                in the grid, an integer
%
%     rmax     - maximum radius that the
%                grid must reach
%
% Outputs:
%
%     phi      - a column vector of polar
%                phi angles, radians
%
%     rmax     - a coluim vector of radii
%
%     L        - sparse Laplacian operator
%                acting on functions defined
%                as vector of values in the
%                same order as the grid
%
% c.brett@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=grid_polar.m>

function [phi,r,L]=grid_polar(ncircles,rmax)

% Check consistency
grumble(ncircles,rmax);

% Get the radius grid
radii=linspace(0,rmax,ncircles);

% Preallocate output arrays
phi=zeros(3*(ncircles-1)^2+2*(ncircles-1)+1,1);
r=zeros(3*(ncircles-1)^2+2*(ncircles-1)+1,1);

% Build grid circles
for n=1:(numel(radii)-1)
    
    % Generate the circle
    phi_current=linspace(0,2*pi,6*n)'; phi_current(end)=[];
    r_current=radii(n+1)*ones(6*n,1);  r_current(end)=[];
    
    % Update the output array
    phi((3*(n-1)^2+2*(n-1)+2):(3*n^2+2*n+1))=phi_current;
    r((3*(n-1)^2+2*(n-1)+2):(3*n^2+2*n+1))=r_current;
             
end

% Build the Laplacian
if nargout>2

    % Make the grid cartesian
    x=r.*cos(phi); y=r.*sin(phi);

    % Run the triangulation
    tri=delaunay(x,y);
    
    % Populate the Laplacian
    L=zeros(numel(r),numel(r));
    for n=1:size(tri,1)
        L(tri(n,1),tri(n,2))=1/norm([x(tri(n,1))-x(tri(n,2)) y(tri(n,1))-y(tri(n,2))],2)^2;
        L(tri(n,1),tri(n,3))=1/norm([x(tri(n,1))-x(tri(n,3)) y(tri(n,1))-y(tri(n,3))],2)^2;
        L(tri(n,2),tri(n,3))=1/norm([x(tri(n,2))-x(tri(n,3)) y(tri(n,2))-y(tri(n,3))],2)^2;
    end

    % Do the miscellaneous cosmetics
    L=L+transpose(L); L=L-diag(sum(L,1)); L=-L/norm(L,2); L=sparse(L);

end
    
end

% Consistency enforcement
function grumble(ncircles,rmax)
if (~isnumeric(ncircles))||(~isreal(ncircles))||...
   (~isscalar(ncircles))||(ncircles<2)||mod(ncircles,1)
    error('ncircles must be a positive real integer greater than one');
end
if (~isnumeric(rmax))||(~isreal(rmax))||...
   (~isscalar(rmax))||(rmax<=0)
    error('rmax must be a positive real number');
end
end

% Буря мглою месит жижу,
% Хочется все время спать:
% Как же сильно ненавижу
% Эту вашу осень, блять.
%
% Yakov Soba4ki

