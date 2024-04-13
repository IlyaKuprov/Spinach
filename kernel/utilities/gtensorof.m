% Returns the g-tensor of the specified spin at the input 
% orientation. Syntax:
%
%            g=gtensorof(spin_system,spin_number)
%
% Parameters:
%
%    spin_number - a positive integer specifying the number
%                  of the spin in the sys.isotopes list
%
% Outputs:
%
%    g - a 3x3 matrix in Bohr magneton units
%
% Note: the same convention (mu=-mu_b*g*S/hbar) is used for
%       the nuclei, meaning that their g-tensors are much 
%       smaller than those of electrons.
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=gtensorof.m>

function g=gtensorof(spin_system,spin_number)

% Check consistency
grumble(spin_system,spin_number);

% Compute the g-tensor
g=-spin_system.inter.zeeman.ddscal{spin_number}*...
   spin_system.inter.gammas(spin_number)*...
   spin_system.tols.hbar/spin_system.tols.muB;

end

% Consistency enforcement
function grumble(spin_system,spin_number)
if ~all(isfield(spin_system,{'inter','tols'}))
    error('spin system object does not contain the required information.');
end
if (~isnumeric(spin_number))||(~isreal(spin_number))||...
   (numel(spin_number)~=1)||(spin_number<1)||(mod(spin_number,1)~=0)
    error('spin_number must be a positive real integer.');
end
if spin_number>spin_system.comp.nspins
    error('spin_number exceeds the number of spins in the system.');
end
end

% Some people's idea of [free speech] is that they are free 
% to say what they like, but if anyone says anything back,
% that is an outrage.
%
% Winston Churchill

