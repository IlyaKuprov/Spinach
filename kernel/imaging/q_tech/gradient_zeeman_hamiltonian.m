% Builds a position-dependent Zeeman Hamiltonian from local fields.
% Syntax:
%
%          H=gradient_zeeman_hamiltonian(spin_system,fields)
%
% Parameters:
%
%  spin_system - Spinach spin system object
%
%       fields - N-by-3 array of local Bx By Bz fields, tesla; one row
%                per spin in spin_system
%
% Outputs:
%
%            H - Hamiltonian contribution, rad/s
%
% Note: the Zeeman tensor stored in spin_system is scaled by the primary
%       magnet field; this function divides by spin_system.inter.magnet
%       to obtain the per-tesla response.
%
% a.arnab@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gradient_zeeman_hamiltonian.m>

function H=gradient_zeeman_hamiltonian(spin_system,fields)

% Check consistency
grumble(spin_system,fields);

% Get matrix dimensions
spin_op=operator(spin_system,'Lz',1);
H=sparse(size(spin_op,1),size(spin_op,2));

% Add local Zeeman terms
for n=1:spin_system.comp.nspins
    zeeman_per_tesla=spin_system.inter.zeeman.matrix{n}/ ...
        spin_system.inter.magnet;
    omega=zeeman_per_tesla*fields(n,:).';
    H=H+omega(1)*operator(spin_system,'Lx',n)+ ...
        omega(2)*operator(spin_system,'Ly',n)+ ...
        omega(3)*operator(spin_system,'Lz',n);
end
end

% Consistency enforcement
function grumble(spin_system,fields)

if ~isfield(spin_system,'inter')||~isfield(spin_system.inter,'magnet')||...
   (spin_system.inter.magnet==0)
    error('spin_system must contain a non-zero primary magnet field.');
end
if (~isnumeric(fields))||(~isreal(fields))||...
   (size(fields,1)~=spin_system.comp.nspins)||(size(fields,2)~=3)||...
   any(~isfinite(fields(:)))
    error('fields must be a finite real N-by-3 array with one row per spin.');
end

end
