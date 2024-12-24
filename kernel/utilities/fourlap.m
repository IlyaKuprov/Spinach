% Returns a Fourier spectral representation of the Laplacian acting
% on a 3D data array. Syntax:
%
%                    L=fourlap(npoints,extents)
%
% Parameters:
% 
%     npoints -  a three-element vector specifying the number of
%                discretization points in each dimension of the
%                3D cube of data that the operator will be acting
%                on, ordered as [X Y Z].
%
%     extents -  a three-element vector specifying axis extents,
%                ordered as [X Y Z].
%
% Outputs:
%
%     L - Fourier spectral Laplacian, a sparse matrix designed to act 
%         on the vectorization of the 3D data array. The dimensions of
%         the data array are assumed to be ordered as [X Y Z].
%
% Note: periodic boundary conditions.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fourlap.m>

function L=fourlap(npoints,extents)

% Check consistency
grumble(npoints,extents);

% Decide the dimensionality
switch numel(npoints)
    
    case 1
        
        % Get differentiation matrices
        [~,Dxx]=fourdif(npoints(1),2);
        
        % Normalize differentiation matrices
        Dxx=(2*pi/extents(1))^2*Dxx;
        
        % Compute the Laplacian
        L=Dxx;
    
    case 2
        
        % Get differentiation matrices
        [~,Dxx]=fourdif(npoints(1),2);
        [~,Dyy]=fourdif(npoints(2),2);
        
        % Normalize differentiation matrices
        Dxx=(2*pi/extents(1))^2*Dxx;
        Dyy=(2*pi/extents(2))^2*Dyy;
        
        % Compute the Laplacian
        L=kron(Dyy,speye(npoints(1)))+...
          kron(speye(npoints(2)),Dxx);
    
    case 3

        % Get differentiation matrices
        [~,Dxx]=fourdif(npoints(1),2);
        [~,Dyy]=fourdif(npoints(2),2);
        [~,Dzz]=fourdif(npoints(3),2);
        
        % Normalize differentiation matrices
        Dxx=(2*pi/extents(1))^2*Dxx;
        Dyy=(2*pi/extents(2))^2*Dyy;
        Dzz=(2*pi/extents(3))^2*Dzz;
        
        % Compute the Laplacian
        L=kron(kron(Dzz,speye(npoints(2))),speye(npoints(1)))+...
          kron(kron(speye(npoints(3)),Dyy),speye(npoints(1)))+...
          kron(kron(speye(npoints(3)),speye(npoints(2))),Dxx);

    otherwise
        
        % Complain and bomb out
        error('incorrect number of spatial dimensions.');
        
end

end

% Consistency enforcement
function grumble(npoints,extents)
if (~isnumeric(npoints))||(~isreal(npoints))||(any(npoints<1))||any(mod(npoints,1)~=0)
    error('npoints must be a three-element vector of positive integers.');
end
if (~isnumeric(extents))||(~isreal(extents))||(any(extents<=0))
    error('extents must be an array of positive real numbers.');
end
end

% Of course it's the same old story. Truth usually
% is the same old story.
%
% Margaret Thatcher

