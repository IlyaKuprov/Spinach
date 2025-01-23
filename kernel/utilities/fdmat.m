% Returns arbitrary-order central finite-difference differentiation 
% matrices (sparse) with unit grid point spacing. Syntax:
%
%                  D=fdmat(dim,nstenc,order,boundary)
%
% Parameters:
%
%     dim       - dimension of the column vector to be
%                 differentiated
%
%     nstenc    - number of points in the finite diffe-
%                 rence stencil
%
%     order     - order of the derivative required
%
%     boundary  - 'wall' fills the edges with sided
%                 finite difference schemes, 'pbc'
%                 assumes periodic boundaries. The
%                 default is 'pbc'.
%
% Outputs:
%
%     D      - finite difference differentiation matrix
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fdmat.m>

function D=fdmat(dim,nstenc,order,boundary)

% Set the default
if ~exist('boundary','var'), boundary='pbc'; end

% Check consistency
grumble(dim,nstenc,order,boundary);

% Preallocate the answer
D=spalloc(dim,dim,dim*nstenc);

% Decide the filling procedure
switch boundary
    
    case 'wall'

        % Edges filled with sided schemes
        for n=1:(nstenc-1)/2
            w=fdweights(n,1:nstenc,order); D(n,1:nstenc)=w(end,:); %#ok<SPRIX>
            D(end-n+1,(end-nstenc+1):end)=((-1)^order)*w(end,end:-1:1); %#ok<SPRIX>
        end
        
        % Middle filled with centered schemes
        stencil=((1-nstenc)/2):((nstenc-1)/2);
        w=fdweights(0,stencil,order);
        for n=((nstenc-1)/2+1):(dim-(nstenc-1)/2)
            D(n,stencil+n)=w(end,:); %#ok<SPRIX>
        end
        
    case 'pbc'
        
        % Wraparound fill with centered schemes
        stencil=((1-nstenc)/2):((nstenc-1)/2);
        w=fdweights(0,stencil,order);
        for n=1:dim
            D(n,mod(stencil+n-1,dim)+1)=w(end,:); %#ok<SPRIX>
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown boundary type.');
        
end

end

% Consistency enforcement
function grumble(dim,nstenc,order,boundary)
if (dim<1)||(nstenc<1)||(order<1)||(mod(dim,1)~=0)||...
   (mod(nstenc,1)~=0)||(mod(order,1)~=0)
    error('all input parameters must be positive integers.');
end
if dim<3
    error('minimum differentiation matrix dimension is 3.');
end
if mod(nstenc,2)~=1
    error('the number of stencil points must be odd.');
end
if order>=nstenc
    error('derivative order must be smaller than the stencil size.');
end
if ~ischar(boundary)
    error('the boundary parameter must be a string.');
end
end

% We may eventually come to realize that chastity 
% is no more a virtue than malnutrition.
%
% Alex Comfort

