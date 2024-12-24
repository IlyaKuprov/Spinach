% Returns the image painted within the Fokker-Planck vector by 
% the user-specified spin state. Syntax:
%
%                   phan=fpl2phan(rho,coil,dims)
%
% Parameters:
%
%     rho    - state vector in Fokker-Planck space
%
%    coil    - observable state vector in Liouville space
%
%    dims    - spatial dimensions of the Fokker-Planck
%              problem, a row vector of integers
%
% Output:
%
%    phan    - the image painted by the specified state
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fpl2phan.m>

function phan=fpl2phan(rho,coil,dims)

% Check consistency
grumble(rho,coil,dims);

% Expose the spin dimension
rho=reshape(rho,[numel(coil) prod(dims)]);

% Compute the observable
phan=coil'*rho;

% Reshape as needed
phan=reshape(phan,dims);

end

% Consistency enforcement
function grumble(rho,coil,dims)
if (~isnumeric(rho))||(size(rho,2)~=1)
    error('rho must be a column vector.');
end
if (~isnumeric(coil))||(size(coil,2)~=1)
    error('coil must be a column vector.');
end
if (~isnumeric(dims))||(size(dims,1)~=1)||...
   (~isreal(dims))||any(dims<1)||any(mod(dims,1)~=0)
    error('dims must be a row vector of real positive integers.');
end
if numel(rho)~=prod(dims)*numel(coil)
    error('numel(rho) does not match space(x)spin Kronecker product.');
end
end

% Malcolm Levitt's most fearsome battle cry, the one that sends the chill
% down the spines of any committee and makes marrow freeze in their bones
% in the expectation of what is to come, is "Erm, excuse me..."

