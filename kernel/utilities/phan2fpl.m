% Projects a spatial intensity distribution into the Fokker-Planck
% space, using it as the image painted by the the spin state supp-
% lied. Syntax:
%
%                      rho=phan2fpl(phan,rho)
%
% Parameters:
%
%        phan  - phantom (the spatial distribution of the
%                amplitude of the specified spin state)
%
%        rho   - Liouville space state vector
%
% Outputs:
%
%        rho   - Fokker-Planck state vector
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=phan2fpl.m>

function rho=phan2fpl(phan,rho)

% Check consistency
grumble(phan,rho);

% Stretch the phantom and kron it with the spin state
rho=kron(phan(:),rho);

end

% Consistency enforcement
function grumble(phan,rho)
if (~isnumeric(rho))||(size(rho,2)~=1)
    error('rho must be a column vector.');
end
if (~isnumeric(phan))||(~ismember(ndims(phan),[1 2 3]))||(~isreal(phan))
    error('the phantom must be a 1D, 2D or 3D array of real numbers.');
end
end

% Q: "How many members of a certain demographic group does
%     it take to perform a specified task?"
%
% A: "A finite number: one to perform the task and the rem-
%     ainder to act in a manner stereotypical of the group
%     in question."

