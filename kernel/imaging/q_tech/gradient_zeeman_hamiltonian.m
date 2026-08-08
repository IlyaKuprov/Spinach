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

% Preallocate the Hamiltonian
H=mprealloc(spin_system,3);

% Add local Zeeman terms
for n=1:spin_system.comp.nspins

    % Get Zeeman response per tesla
    zeeman_per_tesla=spin_system.inter.zeeman.matrix{n}/...
                     spin_system.inter.magnet;
    omega=zeeman_per_tesla*fields(n,:).';

    % Add the Cartesian spin operators
    H=H+omega(1)*operator(spin_system,'Lx',n)+...
        omega(2)*operator(spin_system,'Ly',n)+...
        omega(3)*operator(spin_system,'Lz',n);
end

end

% Consistency enforcement
function grumble(spin_system,fields)

if (~isfield(spin_system,'comp'))||(~isfield(spin_system.comp,'nspins'))||...
   (~isfield(spin_system,'inter'))||(~isfield(spin_system.inter,'magnet'))||...
   (~isfield(spin_system.inter,'zeeman'))||...
   (~isfield(spin_system.inter.zeeman,'matrix'))
    error('spin_system must contain Spinach spin, magnet, and Zeeman information.');
end
if (~isnumeric(spin_system.inter.magnet))||(~isreal(spin_system.inter.magnet))||...
   (~isscalar(spin_system.inter.magnet))||(spin_system.inter.magnet==0)
    error('spin_system must contain a non-zero primary magnet field.');
end
if (~iscell(spin_system.inter.zeeman.matrix))||...
   (numel(spin_system.inter.zeeman.matrix)~=spin_system.comp.nspins)
    error('spin_system Zeeman tensor array must contain one element per spin.');
end
for n=1:spin_system.comp.nspins
    if (~isnumeric(spin_system.inter.zeeman.matrix{n}))||...
       (~isreal(spin_system.inter.zeeman.matrix{n}))||...
       (~isequal(size(spin_system.inter.zeeman.matrix{n}),[3 3]))||...
       any(~isfinite(spin_system.inter.zeeman.matrix{n}(:)))
        error('spin_system Zeeman tensors must be finite real 3x3 matrices.');
    end
end
if (~isnumeric(fields))||(~isreal(fields))||...
   (size(fields,1)~=spin_system.comp.nspins)||(size(fields,2)~=3)||...
   any(~isfinite(fields(:)))
    error('fields must be a finite real N-by-3 array with one row per spin.');
end

end

