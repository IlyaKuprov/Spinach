% Computes axiality and rhombicity of a symmetric 3x3 interaction
% tensor from the corresponding matrix. Syntax:
%
%                     [iso,ax,rh]=mat2axrh(M)
%
% Parameters:
%
%         M   - a real symmetric 3x3 matrix
%
% Outputs:
%
%        iso  - isotropic part of the interaction, defined as
%               (xx+yy+zz)/3 in terms of eigenvaues
%
%         ax  - interaction axiality, defined as zz-(xx+yy)/2
%               in terms of eigenvalues
%
%         rh  - interaction rhombicity, defined as (xx-yy) in
%               terms of eigenvalues
%
% Note: eigenvalues [xx yy zz] are sorted in Mehring order, that
%       is xx<=yy<=zz
%
% Note: Euler angles are not returned because the transformation
%       in question is ill-defined
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mat2axrh.m>

function [iso,ax,rh,eigvals]=mat2axrh(M)

% Check consistency
grumble(M);

% Get the eigenvalues
[~,eigvals]=eig(M); 
eigvals=diag(eigvals);

% Put eigenvalues in Mehring order
eigvals=sort(eigvals,'ascend');

% Get the outputs
iso=mean(eigvals);
ax=eigvals(3)-(eigvals(1)+eigvals(2))/2;
rh=eigvals(1)-eigvals(2);

end

% Consistency enforcement
function grumble(M)
if (~isnumeric(M))||(~isreal(M))||...
   (~all(size(M)==[3 3]))||(~issymmetric(M))
    error('the input must be a real symmetric 3x3 matrix.');
end
end

% People say it would be terrible if we made all 
% girls pretty. I think it would be great.
%
% James Watson

