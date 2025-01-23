% Returns the average of the spin state vector across the 
% spatial dimensions of the sample. Syntax:
%
%                  rho=fpl2rho(rho,dims)
%
% Parameters:
%
%      rho   - Fokker-Planck state vector
%
%     dims   - spatial dimensions of the 
%              Fokker-Planck problem, a
%              vector of positive integers
%
% Outputs:
%
%      rho   - Liouville space state vector
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fpl2rho.m>

function rho=fpl2rho(rho,dims)

% Check consistency
grumble(rho,dims);

% Find out the stack size
stack_size=size(rho,2);

% Expose the spin dimension (no ND sparse support yet)
rho=reshape(full(rho),[size(rho,1)/prod(dims) prod(dims) stack_size]);

% Average over the spatial coordinates
rho=squeeze(sum(rho,2)/prod(dims));

end

% Consistency enforcement
function grumble(rho,dims)
if (~isnumeric(rho))
    error('rho must be a numeric array.');
end
if (~isnumeric(dims))||(size(dims,1)~=1)||...
   (~isreal(dims))||any(dims<1)||any(mod(dims,1)~=0)
    error('dims must be a row vector of real positive integers.');
end
if mod(numel(rho)/prod(dims),1)~=0
    error('numel(rho) does not match space(x)spin Kronecker product.');
end
end

% I would like to die on Mars, just not on impact.
%
% Elon Musk

