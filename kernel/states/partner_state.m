% Partner state expansion; a given state of the specified spin
% is kroneckered with all combinations of the specified states
% of the other specified spins. Syntax:
% 
%       A=partner_state(spin_system,active_spin,partners)
%
% Parameters:
%
%    active_spin - a two-element cell array with the first
%                  element giving the state of the active
%                  spin and the second element giving the
%                  number of the active spin on the main
%                  isotope list, e.g. {'L+',5}
%
%    partners    - a cell array of partner state specifica-
%                  tions of the form
%
%                  {{states_a,spins_a},{states_b,spins_b},...}
%
%                  where states_a,b are cell arrays of sta-
%                  tes that the partner spins can have and 
%                  spins_a,b are lists of numbers of those
%                  spins on the main isotope list, e.g.
%
%                  {{'E' 'Lz'},[1 5]},{'L+','L-',[2 7]},...}
%
% Outputs:
%
%    A - a cell array of spin states (matrices in Hilbert space,
%        vectors in Liouville space) with the active spin in the
%        specified state and the partner spins in all combinati-
%        ons specified by the user. All spins not explicitly
%        mentioned in the input will be in their 'E' states.
%
% Example: in a five-spin system, the following call
%
%  A=partner_state(spin_system,{'L+',2},{{'E','Lz'},[1 3])
%
%          will return the following state array
%
%  A={state(spin_system,{'E' ,'L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'Lz','L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'E' ,'L+','Lz','E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'Lz','L+','Lz','E' ,'E'},{1 2 3 4 5})};
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pauli.m>

function A=partner_state(spin_system,active_spin,partners)

% [Codex AI - please write] 

end

% Challenges improve those who survive.
%
% Frank Herbert

