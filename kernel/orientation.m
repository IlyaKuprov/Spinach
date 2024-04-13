% Anisotropic part of the Hamiltonian for a specific spin system
% orientation. Syntax:
%
%                  H=orientation(Q,euler_angles)
% 
% Arguments:
%
%   Q             -  rotational basis as returned by 
%                    hamiltonian.m function
%
%   euler_angles  -  a 1x3 vector specifying Euler 
%                    angles (radians) relative to the
%                    input orientation
%
% Output:
%
%     H    -    anisotropic part of the Hamiltonian 
%               for the specified Euler angles
%
% Note: this function may be used in both Hilbert and Liouville 
%       space because the H -> [H, ] adjoint map is linear.
%
% TODO: efficient sparse summation.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=orientation.m>

function H=orientation(Q,euler_angles)

% Check consistency
grumble(Q,euler_angles);

% Get started
H=sparse(0);

% Loop over the spherical ranks
for r=1:numel(Q)
    
    % Compute the Wigner matrix
    D=wigner(r,euler_angles(1),...
               euler_angles(2),...
               euler_angles(3));

    % Update the Hamiltonian
    for k=1:(2*r+1)
        for m=1:(2*r+1)
            if nnz(Q{r}{k,m})>0
                H=H+D(k,m)*Q{r}{k,m};
            end
        end
    end
    
end
    
% Clean up the Hamiltonian
H=(H+H')/2;    

end

% Consistency enforcement
function grumble(Q,euler_angles)
if ~iscell(Q)
    error('Q parameter must be a cell array.');
end
if (~isnumeric(euler_angles))||...
   (~isreal(euler_angles))||...
   (numel(euler_angles)~=3)
    error('euler_angles must be a three-element real vector.')
end
end

% "Can you explain how Spinach calculates protein structures?"
% "But surely we can already compute NMR spectra with DFT?"
% "Can you say what is actually novel in any of this?"
%
% Ulrike Salzner's questions to IK at the 2015
% ERC Panel. The proposal was not funded, stal-
% ling Spinach development for at least a year.

