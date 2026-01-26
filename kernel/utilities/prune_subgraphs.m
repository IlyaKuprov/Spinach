% Removes subgraphs that are contained entirely within
% other subgraphs. Syntax:
%
%         subgraphs=prune_subgraphs(subgraphs)
%
% Parameters:
%
%    subgraphs - [ngraphs x nspins] logical array
%                with 1 when a spin belongs to a
%                subgraph and 0 otherwise
%
% Outputs:
%
%    subgraphs - [ngraphs x nspins] logical array
%                with 1 when a spin belongs to a
%                subgraph and 0 otherwise
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=prune_subgraphs.m>

function subgraphs=prune_subgraphs(subgraphs)

% Check consistency
grumble(subgraphs);

% Ignore trivial cases
if (size(subgraphs,1)<2)||...
   (size(subgraphs,2)<2)
    return;
end

% Count spins in each subgraphs
spin_counts=sum(subgraphs,2);

% Get subgraph overlap matrix
C=subgraphs*transpose(subgraphs);

% Check for supersets    
supersets=((C==spin_counts)&...
           (transpose(spin_counts)>spin_counts));

% Do the pruning
subgraphs=subgraphs(~any(supersets,2),:);

end

% Consistency enforcement
function grumble(subgraphs)
if ~islogical(subgraphs)
    error('subgraphs must be a logical array.');
end
end

% Evil people always support each other;
% that is their chief strength.
%
% Alexander Solzhenytsyn

