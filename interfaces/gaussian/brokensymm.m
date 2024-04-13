% Exchange coupling estimation from a pair of DFT logs using 
% Yamaguchi equation. The notation is:
%
%                        H=-2J*(Sa.Sb)
%
% Syntax:
%
%             J=brokensymm(props_sing,props_trip)
%
% Parameters:
%
%     props_sing  - the output of gparse for the singlet 
%                   state of the biradical
%
%     props_trip  - the output of gparse for the triplet 
%                   state of the biradical
%
% Outputs:
%
%     J           - an order-of-magnitude (really rough)
%                   estimate of exchange coupling, Hz
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=brokensymm.m>

function J=brokensymm(props_sing,props_trip)

% Check consistency
grumble(props_sing,props_trip);

% Eq 6 in https://doi.org/10.1063/1.5144696
J=(props_trip.energy-props_sing.energy)/...(
  (props_sing.s_sq-props_trip.s_sq);

% Convert from Hartree to Hz
J=6.57968974479e15*J;

end

% Consistency enforcement
function grumble(props_sing,props_trip)
if ~isfield(props_sing,'energy')
    error('props_sing.energy field is missing');
end
if ~isfield(props_trip,'energy')
    error('props_trip.energy field is missing');
end
if ~isfield(props_sing,'s_sq')
    error('props_sing.s_sq field is missing');
end
if ~isfield(props_trip,'s_sq')
    error('props_trip.s_sq field is missing');
end
end

% "There is nothing new to be discovered in physics now. All
%  that remains is more and more precise measurement."
% 
% Lord Kelvin, in 1900, addressing the British
% Association for the Advancement of Science

