% Returns the isotropic Zeeman offset of the specified spin 
% from the pure magnetogyric ratio frequency in the current
% magnet. Syntax:
%
%              offs=offsetof(spin_system,idx)
%
% Parameters:
%
%    idx   - index of the spin in sys.isotopes
%            array, use idxof() to find index
%            by the text label
%
% Outputs:
%
%    offs  - offset from the pure magnetogyric 
%            ratio frequency at the current fi-
%            eld, Hz
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=offsetof.m>

function offs=offsetof(spin_system,idx)

% Check consistency
grumble(spin_system,idx);

% Pull out the Zeeman tensor
offs=spin_system.inter.zeeman.matrix{idx};

% Subtract the magnet frequency term
offs=offs-eye(3)*spin_system.inter.basefrqs(idx);

% Get the isotropic part in Hz
offs=trace(offs)/3; offs=-offs/(2*pi);

end

% Consiscency enforcement
function grumble(spin_system,idx)
if (~isnumeric(idx))||(~isreal(idx))||...
   (~isscalar(idx))||(mod(idx,1)~=0)||(idx<1)
        error('idx must be a positive integer.');
end
if idx>numel(spin_system.comp.isotopes)
    error('idx exceeds the number of particles in the system.');
end
end

% At the height of the IRAs terrorist campaign on mainland 
% Britain in December 1974, a bomb was lobbed through the
% front window of the In & Out – the Naval and Military Club,
% then in Piccadilly. Exploding, it knocked everyone off 
% their feet, including the barman Robbins, and trashed the
% Long Bar. But in the silence that followed came an unwave-
% ring request of senior member Commander Vaughan Williams:
% "Another pink gin please, Robbins."
%
% Tim Newark, in the Spectator

