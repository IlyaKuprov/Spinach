% Generates sets of vector-covector pairs for the parallel
% implementation of the time propagation algorithm described in
% http://dx.doi.org/10.1063/1.3679656 (Equation 9). Syntax:
%
%         [vec,cov]=svd_shrink(spin_system,rho,tol)
%
% Parameters:
%
%         rho   -   density matrix
%
%         tol   -   singluar value drop tolerance
%
% Outputs:
%
%         vec   -   vectors as columns of a matrix
%
%         cov   -   covectors as columns of a matrix
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de 
%
% <https://spindynamics.org/wiki/index.php?title=svd_shrink.m>

function [vec,cov]=svd_shrink(spin_system,rho,tol)

% Check consistency
grumble(rho);

% Run the singular value decomposition
[vec,S,cov]=svd(full(rho)); S=diag(S);

% Get the drop mask
drop_mask=(S<tol);

% Update the user
report(spin_system,['dropped ' num2str(nnz(drop_mask)) ' insignificant '...
                    'vector-covector pairs from the density matrix.']);

% Eliminate small singular values
vec(:,drop_mask)=[]; cov(:,drop_mask)=[]; S(drop_mask)=[];

% Spread the coefficients
vec=vec*diag(sqrt(S));
cov=cov*diag(sqrt(S));

end

% Consistency enforcement
function grumble(rho)
if (~isnumeric(rho))||(size(rho,1)~=size(rho,2))
    error('rho must be a square matrix.');
end
end

% If everyone likes your research, you can be certain that
% you have not done anything important. That is the first
% thing to grasp. Conflict goes with the territory.
%
% Andrew Oswald

