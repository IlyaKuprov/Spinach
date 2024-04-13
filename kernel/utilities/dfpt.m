% Graph partitioning module. Analyzes the system connectivity graph and
% creates a list of all connected subgraphs of up to the user-specified
% size by crawling the graph in all available directions. Syntax:
%
%                 subgraphs=dfpt(conmatrix,max_sg_size)
%
% Parameters:
%
%         conmatrix   - [nspins x nspins] matrix with 1 for 
%                       connected spins and 0 elsewhere
%
%       max_sg_size   - maximum connected subgraph size
%
% Outputs: 
%
%         subgraphs   - [n_subgraphs x nspins] matrix; each 
%                       row contains 1 for spins that belong
%                       to the subgraph and 0 for spins that
%                       do not
%
% i.kuprov@soton.ac.uk
% stst@uw.edu
%
% <https://spindynamics.org/wiki/index.php?title=dfpt.m>

function subgraphs=dfpt(conmatrix,max_sg_size)

% Check consistency
grumble(conmatrix,max_sg_size);

% Start at each spin
nspins=size(conmatrix,2);
subgraphs=uint32(1:nspins)';

% Crawl the graph
for sg_size=2:max_sg_size
    
    % Preallocate grown set
    grown_set=cell(size(subgraphs,1),1);
    
    % Loop over subgraphs
    for n=1:size(subgraphs,1)
        
        % Get current subgraph
        subgraph=subgraphs(n,:);
        
        % Find spins reachable from current subgraph
        neighbours=any(conmatrix(:,subgraph),2);
        neighbours(subgraph)=false();
        neighbours=uint32(find(neighbours));
        
        % Grow in every direction
        if isempty(neighbours)
            
            % Isolated subgraphs get a dummy index
            grown_set{n}=[subgraph subgraph(end)];
            
        else
            
            % Subgraphs with neighbours get neighbours
            grown_set{n}=[repmat(subgraph,[numel(neighbours) 1]) neighbours];
            
        end
        
    end
    
    % Merge grown set
    subgraphs=cell2mat(grown_set);
    subgraphs=sort(subgraphs,2);
    subgraphs=unique(subgraphs,'rows');
        
end

% Get sparse row index
row_index=uint32(1:size(subgraphs,1))';
row_index=repmat(row_index,[1 max_sg_size]);
row_index=row_index(:);

% Get sparse column index
nsg=size(subgraphs,1); subgraphs=subgraphs(:);

% Return subgraph array as a sparse logical matrix
subgraphs=sparse(row_index,subgraphs,1,nsg,nspins);
subgraphs=logical(subgraphs);
        
end

% Consistency enforcement
function grumble(conmatrix,max_sg_size)
if ~islogical(conmatrix)
    error('conmatrix must be a logical matrix.');
end
if size(conmatrix,1)~=size(conmatrix,2)
    error('conmatrix must be square.');
end
if (~isscalar(max_sg_size))||(~isnumeric(max_sg_size))||...
   (~isreal(max_sg_size))||(max_sg_size<1)||...
   (mod(max_sg_size,1)~=0)
    error('max_sg_size must be a real positive integer.');
end
end

% Psychologists call this "cognitive dissonance" - the ability to make a
% compelling, heartfelt case for one thing while doing another. Being able
% to pull off this sort of trick is an essential skill in many professions:
% "even if his message bore no relation to his actions, it expressed pre-
% cisely and succinctly what he should have been doing".
%
% The Economist

