% Chemical reaction generator builder. Syntax:
%
%              G=react_gen(spin_system,reaction)
%
% Parameters:
%
%    reaction.reactants - a vector of integers specifying
%                         which parts declared in the in-
%                         put (chem.parts) are reactants
%
%    reaction.products  - a vector of integers specifying
%                         which parts declared in the in-
%                         put (chem.parts) are products
%
%    reaction.matching  - a matrix with two columns, spe-
%                         cifying which spin in the reac-
%                         tants list (left column) becom-
%                         es which spin in the product 
%                         list (right column)
%
% Outputs:
%
%    G - a matrix mapping each state of the
%        reactant state space into its desti-
%        nation in the product state space;
%        a cell array, one per reactant
%
% Notes: pilot implementation with slow indexing, needs
%        a high-performance overhaul
%
% i.kuproprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=react_gen.m>

function [GD,GF]=react_gen(spin_system,reaction)

% Check consistency
grumble(spin_system,reaction);

% Inform the user
report(spin_system,'building reaction generators...');

% Get drain and fill generators started
GD=cell([numel(reaction.reactants) 1]);
for n=1:numel(GD), GD{n}=zeros([0 3]); end
GF=cell([numel(reaction.reactants) 1]);
for n=1:numel(GF), GF{n}=zeros([0 3]); end

% Loop over the basis set
for n=1:size(spin_system.bas.basis,1)
    
    % Extract the state
    source_state=spin_system.bas.basis(n,:);

    % Find participating spins
    [~,spins_involved]=find(source_state);

    % Build reaction generators
    if ~isempty(spins_involved)

        % Determine the host substance
        host_subst=cellfun(@(x)all(ismember(spins_involved,x)),...
                           spin_system.chem.parts);
        [~,host_subst]=find(host_subst);

        % Double-check indexing
        if numel(host_subst)~=1
            error('basis set indexing problem.');
        end

        % Determine host substance type
        this_is_reactants=ismember(host_subst,reaction.reactants);
        this_is_products=ismember(host_subst,reaction.products);

        % Double-check indexing
        if this_is_reactants&&this_is_products
            error('reactant/product indexing problem.');
        end

        % Only reactants react
        if this_is_reactants

            % Add to reactant drain generator
            idx=find(reaction.reactants==host_subst);
            GD{idx}=[GD{idx}; [n n -1]];
        
            % Build the destination state
            destin_state=zeros(1,numel(source_state));
            destin_state(reaction.matching(:,2))=source_state(reaction.matching(:,1));
            destin_state=sparse(destin_state);

            % Look for the destination state in the basis set and double-check indexing
            [destin_exists,destin_index]=ismember(destin_state,spin_system.bas.basis,'rows');
            if numel(destin_index)>1, error('invalid basis set specification'); end

            % Build product fill generator
            if destin_exists

                % Find participating spins
                [~,spins_involved]=find(destin_state);

                % Determine the host substance
                host_subst=cellfun(@(x)all(ismember(spins_involved,x)),...
                                   spin_system.chem.parts);
                [~,host_subst]=find(host_subst);

                % Double-check indexing
                if numel(host_subst)~=1
                    error('basis set indexing problem.');
                end
                
                % Add to product fill generator
                GF{idx}=[GF{idx}; [destin_index n 1]];

            end

        end

    end

end

% Convert to sparse matrices
dim=size(spin_system.bas.basis,1);
for n=1:numel(GD)
    GD{n}=sparse(GD{n}(:,1),...
                 GD{n}(:,2),...
                 GD{n}(:,3),dim,dim);
end
for n=1:numel(GF)
    GF{n}=sparse(GF{n}(:,1),...
                 GF{n}(:,2),...
                 GF{n}(:,3),dim,dim);
end

end

% Consistency enforcement
function grumble(spin_system,reaction)
if ~isempty(intersect(reaction.reactants,...
                      reaction.products))
    error('reactants and products must contain different substances.');
end
if numel(cell2mat(spin_system.chem.parts(reaction.reactants)))~=...
   numel(cell2mat(spin_system.chem.parts(reaction.products)))
    error('number of spins not the same either side of the reaction arrow.');
end
if (~isempty(setdiff(reaction.matching(:,1)',...
                     cell2mat(spin_system.chem.parts(reaction.reactants)))))||...
   (~isempty(setdiff(reaction.matching(:,2)',...
                     cell2mat(spin_system.chem.parts(reaction.products)))))
    error('matching map and part specification are not consistent.');
end
end

% "Once, years ago, making small talk with Elizabeth II, I asked
%  her if it was true that many peers attending her coronation in
%  1953 had taken sandwiches into Westminster Abbey hidden inside
%  their coronets. 'Oh, yes' - she said. 'They were in the Abbey
%  for something like six hours, you know. The Archbishop of Can-
%  terbury even had a flask of brandy tucked inside his cassock.'
%  Apparently, His Grace offered Her Majesty a discreet nip, but
%  she declined."
%
% Gyles Brandreth

