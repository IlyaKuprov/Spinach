% Extracts the 3x3 coupling tensor between a pair of spins back 
% from the spin_system data structure. Syntax:
%
%                A=get_coupling(spin_system,n,k)
%
% Parameters:
%
%   n,k  - indices of the two spins as they appear 
%          in spin_system.comp.isotopes
%
% Outputs:
%
%     A  - 3x3 coupling tensor in rad/s
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=get_coupling.m>

function A=get_coupling(spin_system,n,k)

% Check consistency
grumble(spin_system,n,k);

% Pull forward and backward coupling
A=spin_system.inter.coupling.matrix{n,k};
B=spin_system.inter.coupling.matrix{k,n};

% Fill in empties
if isempty(A), A=zeros([3 3]); end
if isempty(B), B=zeros([3 3]); end

% Add up
A=A+B;

end

% Consistency enforcement
function grumble(spin_system,n,k)
if (~isfield(spin_system,'inter'))||...
   (~isfield(spin_system.inter,'coupling'))
    error('spin_system object does not contain coupling information.');
end
if (~isnumeric(n))||(~isscalar(n))||(~isreal(n))||(n<1)||(mod(n,1)~=0)||...
   (~isnumeric(k))||(~isscalar(k))||(~isreal(k))||(k<1)||(mod(k,1)~=0)
    error('n and k must be positive real integers.');
end
end

% Единственное, что я понимаю в арбузах - это если я по 
% арбузу постучала, а из него постучали в ответ, то вот
% это прям однозначно плохой арбуз.
%
% Russian internet folklore

