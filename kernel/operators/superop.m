% Product, commutation, and anticommutation superoperators in the 
% spherical tensor basis set for the specified substance. Syntax:
%
%              A=p_superop(spin_system,opspec,side)
%
% Arguments:
%
%     opspec - Spinach operator specification described in Sec-
%              tions 2.1 and 3.3 of the following paper:
%
%                http://dx.doi.org/10.1016/j.jmr.2010.11.008
%
%     side   - 'left' or 'right' causes the function to return
%              a product superoperator corresponding to a pro-
%              duct from that side; 'comm' or 'acomm' results
%              in commutation and anticommutation superopera-
%              tor respectively.
%
% Outputs:
%
%          A - a three-column array of row indices (first col-
%              umn), column indices (second column) and values
%              (third column)
%
% Note: direct calls to this function are not usually required,
%       use the friendlier operator() function.
%
% Note: the superoperator is returned in XYZ sparse format.
%
% ilya.kuprov@weizmann.ac.il
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=p_superop.m>

function A=p_superop(spin_system,opspec,side)

% Recursions
switch side

    case 'comm'

        % Efficient left and right commutator parts,
        % skipping elements that would cancel anyway
        A=p_superop(spin_system,opspec, 'leftofcomm');
        B=p_superop(spin_system,opspec,'rightofcomm');

        % XYZ sparse indexed subtraction
        B(:,3)=-B(:,3); A=[A; B]; return;

    case 'acomm'

        % Full left and right anticommutator parts
        A=p_superop(spin_system,opspec, 'left');
        B=p_superop(spin_system,opspec,'right');

        % XYZ sparse indexed addition
        A=[A; B]; return;

    case {'left','right','leftofcomm','rightofcomm'}

        % Validate the input
        grumble(spin_system,opspec);

    otherwise

        % Complain and bomb out
        error('unrecognised side specification.');

end

% Index relevant spins
opspec_mask=logical(opspec);
active_spins=find(opspec_mask);

% Find out which substance we are in
subst=which_subst(spin_system,active_spins);

% Get relevant spins from substance list
spins_in_subst=spin_system.chem.parts{subst};

% Get the relevant multiplicities from the spin list
mults_in_subst=spin_system.comp.mults(spins_in_subst);

% Start state tables and the
% structure coefficient table
source=cell(0,numel(active_spins));
destin=cell(0,numel(active_spins));
struct=cell(0,numel(active_spins));

% Get action tables
for n=active_spins

    % Spin multiplicity
    mult=mults_in_subst(n);

    % 1-base indexing
    table_idx=opspec(n)+1;
   
    % Extract pages
    switch side

        case {'left','leftofcomm'}
            
            % Extract left product pages from Lie structure tables
            pt=squeeze(spin_system.bas.lpst{mult}(table_idx,:,:));

        case {'right','rightofcomm'}

            % Extract right product pages from Lie structure tables
            pt=squeeze(spin_system.bas.rpst{mult}(table_idx,:,:));
        
        otherwise

            % Complain and bomb out
            error('invalid side specification.');

    end
    
    % Convert product action table to indices
    [destin{n},source{n},struct{n}]=find(pt);
    
    % Spinach uses 0 index for the unit matrix
    source{n}=source{n}-1; destin{n}=destin{n}-1;
    
end
    
% Get structure coeffs for the sub-algebra
from=source{1}; to=destin{1}; coeff=struct{1};
for n=2:numel(active_spins)
    from=[repelem(from,height(source{n}),1) repmat(source{n},height(from),1)];
    to=  [repelem(to,height(destin{n}),1)   repmat(destin{n},height(to),1)  ];
    coeff=kron(coeff,struct{n});
end

% Get basis columns corresponding to the relevant spins
basis_cols=spin_system.bas.basis{subst}(:,opspec_mask);

% For commutation superoperators, remove commuting paths
if ismember(side,{'leftofcomm','rightofcomm'})
    idx=(sum(from,2)==0)|(sum(to,2)==0);
    from(idx,:)=[]; to(idx,:)=[]; coeff(idx,:)=[];
end

% Get the answer going
A=cell(height(from),1);

% Loop over source states
for n=1:size(from,1)
    
    % Retrieve the source subspace
    source_subsp_idx=true(height(basis_cols),1);
    for m=1:width(from)
        source_subsp_idx=and(source_subsp_idx,...
                             basis_cols(:,m)==from(n,m));
    end
    source_subsp=spin_system.bas.basis{subst}(source_subsp_idx,...
                                              ~opspec_mask);
    source_subsp_idx=find(source_subsp_idx);
    
    % Get source subspace dimension
    subsp_dim=size(source_subsp,1);
    
    % Skip empty subspaces
    if subsp_dim>0
    
        % Retrieve the destination subspace
        destin_subsp_idx=true(height(basis_cols),1);
        for m=1:width(to)
            destin_subsp_idx=and(destin_subsp_idx,...
                                 basis_cols(:,m)==to(n,m));
        end
        destin_subsp=spin_system.bas.basis{subst}(destin_subsp_idx,...
                                                  ~opspec_mask);
        destin_subsp_idx=find(destin_subsp_idx);
        
        % Fill the operator
        if isequal(source_subsp,destin_subsp)
            
            % If the subspaces fully match, use the raw indices
            A{n}=[source_subsp_idx destin_subsp_idx repelem(coeff(n),subsp_dim,1)];
            
        else
        
            % Otherwise, use brute-force state-by-state matching
            [does_it_go_anywhere,where_it_goes_if_it_does]=ismember(source_subsp, ...
                                                                    destin_subsp,'rows');
            A{n}=[source_subsp_idx(does_it_go_anywhere)                           ...
                  destin_subsp_idx(where_it_goes_if_it_does(does_it_go_anywhere)) ...
                  repelem(coeff(n),nnz(does_it_go_anywhere),1)];
        
        end
        
    end
    
end

% Remove empty cells
A(cellfun(@isempty,A))=[];

% Merge cells
if isempty(A)
    A=[1 1 0];
else
    A=cell2mat(A);
end

% Get the global multi-substance basis index offset
idx_offset=sum(spin_system.bas.nstates(1:(subst-1)));

% Apply the offset
A(:,1)=A(:,1)+idx_offset;
A(:,2)=A(:,2)+idx_offset;

end

% Consistency enforcement
function grumble(spin_system,opspec)

% Basis information must exist
if ~isfield(spin_system,'bas')
    error('basis set information is missing, call basis() first.');
end

% Spherical tensors in Liouville space only
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('p_superop() requires sphten-liouv formalism.');
end

% Opspec must be sensible and physically valid
if (~isnumeric(opspec))||(~isrow(opspec))||any(mod(opspec,1)~=0)
    error('opspec must be a row vector of integers.');
end
if numel(opspec)~=spin_system.comp.nspins
    error('numel(opspec) must be equal to the number of spins.');
end
if any((opspec+1)>spin_system.comp.mults.^2)
    error('physically impossible state requested in opspec.');
end
if nnz(opspec)==0
    error('for all-zero opspec call unit_state() instead.'); 
end

end

% My philosophy, in essence, is the concept of man as a heroic being,
% with his own happiness as the moral purpose of his life, with pro-
% ductive achievement as his noblest activity, and reason as his only
% absolute.
%
% Ayn Rand

