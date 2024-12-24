% Returns the chemical shifts of every spin in the system relative
% to the carrier frequency in the current magnet. Syntax:
%
%               [cs_ppm,cs_hz]=chemshifts(spin_system)
%
% Outputs:
%
%    cs_ppm  - chemical shifts in ppm
%
%    cs_hz   - chemical shifts in Hz
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=chemshifts.m>

function [cs_ppm,cs_hz]=chemshifts(spin_system)

% Check consistency
grumble(spin_system);

% Preallocate outputs
cs_ppm=zeros(spin_system.comp.nspins,1);
cs_hz=zeros(spin_system.comp.nspins,1);

% Fill the outputs
for n=1:spin_system.comp.nspins
    
    % Get isotropic Zeeman frequencies
    iso=trace(spin_system.inter.zeeman.matrix{n})/3;
    
    % Subtract carrier
    iso=iso-spin_system.inter.basefrqs(n);
    
    % Convert into ppm
    cs_ppm(n)=1e6*iso/spin_system.inter.basefrqs(n);
    
    % Convert into Hz
    cs_hz(n)=-iso/(2*pi);
    
end

end

% Consistency enforcement
function grumble(spin_system)
if (~isfield(spin_system,'comp'))||(~isfield(spin_system,'inter'))
    error('the spin system object does not contain the required information.');
end
end

% I have never understood why it is 'greed' to want to 
% keep the money you have earned, but not greed to want
% to take somebody else's money.
%
% Thomas Sowell

