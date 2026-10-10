% Compile product-row maps for an explicit chemical reaction. Syntax:
%
%                   maps=react_gen(spin_system,reaction)
%
% Parameters:
%
%    spin_system - Spinach object with a compiled sphten-liouv basis
%
%    reaction    - validated chem.reactions record from create()
%
% Outputs:
%
%    maps        - cell array, one per product occurrence; each matrix
%                  has global destination indices in its first column
%                  and one global source index per reactant occurrence
%                  in the remaining columns
%
% Product-only spins must be at identity. Unmatched source spins are
% traced out. Missing source descriptors are counted and reported;
% no Cartesian product of the reactant bases is materialised. Repeated
% spin-bearing reactants or products with matched spins require occurrence-resolved
% matching and are rejected by this two-column matching interface.
%
% i.kuproprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=react_gen.m>

function maps=react_gen(spin_system,reaction)

% Check consistency
grumble(spin_system,reaction);

% Enumerate each product block once per stoichiometric occurrence
maps=cell(size(reaction.products));
for n=1:numel(reaction.products)

    % Require unpolarised arrival on spins without a source
    product=reaction.products(n);
    spins=spin_system.chem.parts{product};
    descr=spin_system.bas.basis{product};
    present=~any(descr(:,~ismember(spins,reaction.matching(:,2))),2);
    sources=zeros(size(descr,1),numel(reaction.reactants));
    missing=false(size(present));
    for k=1:numel(reaction.reactants)

        % Pull the product descriptor back through the atom matching
        reactant=reaction.reactants(k);
        source_spins=spin_system.chem.parts{reactant};
        [on_source,source_col]=ismember(reaction.matching(:,1),source_spins);
        [on_product,product_col]=ismember(reaction.matching(:,2),spins);
        matched=on_source&on_product;
        [rows,cols,values]=find(descr(:,product_col(matched)));
        source_col=source_col(matched);
        source_descr=sparse(rows,source_col(cols),values,size(descr,1),numel(source_spins));

        % Locate source rows without constructing a reactant tensor basis
        if isempty(source_spins)
            found=true(size(present)); index=ones(size(present));
        else
            [found,index]=ismember(source_descr,spin_system.bas.basis{reactant},'rows');
        end
        missing=missing|~found;
        sources(:,k)=spin_system.bas.offsets(reactant)+index;
    end

    % Report truncation separately from unmatched product-spin orders
    report(spin_system,['reaction product ' num2str(product) ': ' ...
           num2str(nnz(present&missing)) ' product rows without a source']);
    present=present&~missing;
    maps{n}=[spin_system.bas.offsets(product)+find(present) sources(present,:)];
end

end

% Consistency enforcement
function grumble(spin_system,reaction)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('Spinach:react_gen:formalism','reaction maps currently require sphten-liouv formalism.');
end
if ~isstruct(reaction)||~isscalar(reaction)||...
   ~all(isfield(reaction,{'reactants','products','matching'}))
    error('Spinach:react_gen:record','reaction must be a validated chem.reactions record.');
end
for n=unique(reaction.reactants)
    if nnz(reaction.reactants==n)>1&&...
       any(ismember(reaction.matching(:,1),spin_system.chem.parts{n}))
        error('Spinach:react_gen:repeatedMatching',...
              'matched repeated reactants require occurrence-resolved matching, not a two-column spin map.');
    end
end
for n=unique(reaction.products)
    if nnz(reaction.products==n)>1&&...
       any(ismember(reaction.matching(:,2),spin_system.chem.parts{n}))
        error('Spinach:react_gen:repeatedProductMatching',...
              'matched repeated products require occurrence-resolved matching, not a two-column spin map.');
    end
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

