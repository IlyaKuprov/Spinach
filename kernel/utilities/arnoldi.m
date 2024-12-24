% Arnoldi procedure for the creation of an orthonormal Krylov basis
% from repeated action by an operator on a vector. The procedure is
% numerically unstable and must be used with caution. Syntax:
%
%                     [V,H]=arnoldi(Op,v0,niter)
%
% Parameters:
%
%    Op     - function handle taking in a column vector 
%             and returning another column vector
%
%    v0     - starting vector of the Arnoldi process
%
%    nsteps - number of iterations to take; the Krylov
%             subspace will be nsteps+1 dimensional
%
% Outputs:
%
%    V - a matrix containing the orthonormal basis vec-
%        tors of the Krylov subspace in columns
%
%    H - extended Hessenberg matrix
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=arnoldi.m>

function [V,H]=arnoldi(Op,v0,niter)

% Check consistency
grumble(Op,v0,niter);

% Preallocate outputs
V=zeros(numel(v0),niter+1,'like',1i);
H=zeros(niter+1,niter,'like',1i);

% First vector
V(:,1)=v0/norm(v0,2);

% Iteration loop
for n=1:niter
    
    % Get the next vector
    V(:,n+1)=Op(V(:,n));
    
    % Gram-Schmidt loop
    for k=1:n
        
        % Get the scalar product
        H(k,n)=V(:,k)'*V(:,n+1);
        
        % Remove its contribution
        V(:,n+1)=V(:,n+1)-H(k,n)*V(:,k);
        
    end
    
    % Compute the norm
    H(n+1,n)=norm(V(:,n+1),2);
    
    % Divide out the norm
    V(:,n+1)=V(:,n+1)/H(n+1,n);
    
end

end

% Consistency enforcement
function grumble(Op,v0,niter)
if ~isa(Op,'function_handle')
    errror('Op must be a function handle.');
end
if (~isnumeric(v0))||(~iscolumn(v0))
    error('v0 must be a column vector');
end
if (~isnumeric(niter))||(~isreal(niter))||...
   (~isscalar(niter))||(mod(niter,1)~=0)||niter<0
    error('niter must be a non-negative real integer');
end
end

% Life cannot just be about solving one sad problem after 
% another. There need to be things that inspire you, that
% make you glad to wake up in the morning and be part of 
% humanity. That is why we did it.
%
% Elon Musk, about launching his car 
% in the general direction of Mars.

