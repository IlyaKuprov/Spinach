% Adiabatic elimination in Liouville space, implements 
% Section 6.1 of Kuprov's book. Syntax:
%
%    [L,R]=adelim(spin_system,L,fast_idx,slow_idx)
%
% Parameters:
%
%   L        - Liouvillian in sphten-liouv formalism,
%              fast subbsystem must be dissipative
%
%   fast_idx - a vector of integers specifying which
%              states in the basis involve the fast
%              subsystem in any way
%
%   slow_idx - a vector of integers specifying which
%              states in the basis only involve the
%              slow subsystem
%
% Outputs:
%
%   L        - projection of the original Liouvillian
%              into the slow subspace, inheriting any
%              coherent and dissipative dynamics that
%              the user previously had there
%
%   R        - the extra relaxation superoperator on-
%              ce the fast subspace is adiabatically
%              eliminated
%
% Note: the function needs sphten-liouv formalism because
%       there the basis states are attributable to indivi-
%       dual spins.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=adelim.m>

function [L,R]=adelim(spin_system,L,fast_idx,slow_idx)

% Check consistency
grumble(spin_system,L,fast_idx,slow_idx);

% Get unit operator
U=unit_oper(spin_system);

% Get projectors (faster than indexing)
P_slow=U(slow_idx,:); P_fast=U(fast_idx,:);

% Get projections (faster than indexing)
L01=P_slow*L*P_fast'; L10=P_fast*L*P_slow'; 
L11=P_fast*L*P_fast'; L00=P_slow*L*P_slow';

% Eq 6.2 in Kuprov's book
R=1i*L01*(L11\L10); L=L00;

end

% Consistency enforcement
function grumble(spin_system,L,fast_idx,slow_idx)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('adiabatic elimination is only available for sphten-liouv formalism.');
end
if (~isnumeric(L))||(size(L,1)~=size(L,2))
    error('L must be a square matrix.');
end
if (~isnumeric(fast_idx))||(~isvector(fast_idx))||...
   (~isnumeric(slow_idx))||(~isvector(slow_idx))
    error('fast_idx and slow_idx must be vectors.');
end
if ~isempty(intersect(fast_idx,slow_idx))
    error('fast_idx and slow_idx cannot have common elements.');
end
if (numel(fast_idx)+numel(slow_idx))~=size(L,1)
    error('dimensions or indices internally inconsistent.');
end
end

% I have a foreboding of an America in my children's or grand-
% children's time, when the United States is a service and in-
% formation economy; when nearly all the manufacturing indust-
% ries have slipped away to other countries; when awesome tech-
% nological powers are in the hands of very few, and no one
% representing the public interest can even grasp the issues;
% when the people have lost the ability to set their own agen-
% das or knowledgeably question those in authority; when, clut-
% ching our crystals and nervously consulting our horoscopes,
% our critical faculties in decline, unable to distinguish bet-
% ween what feels good and what's true, we slide, almost with-
% out noticing, back into superstition and darkness. The dumb-
% ing down of America is most evident in the slow decay of sub-
% stantive content in the enormously influential media, the 30
% second sound bites (now down to 10 seconds or less), lowest
% common denominator programming, credulous presentations on
% pseudoscience and superstition, but especially a kind of ce-
% lebration of ignorance.
%
% Carl Sagan, 1995

