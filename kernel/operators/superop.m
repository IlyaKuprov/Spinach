% Sided product superoperator in the spherical tensor basis set. Returns
% superoperators corresponding to right or left multiplication of a den-
% sity matrix by a user-specified operator. Syntax:
%
%                  A=superop(spin_system,opspec,side)
%
% Arguments:
%
%     opspec - Spinach operator specification described in Sections 2.1
%              and 3.3 of the following paper:
%
%                     http://dx.doi.org/10.1016/j.jmr.2010.11.008
%
%     side   - 'left' or 'right' causes the function to return a product
%              superoperator corresponding to a product from that side;
%              'comm' or 'acomm' results in commutation and anticommuta-
%              tion superoperator respectively.
%
% Outputs:
%
%          A - a three-column array of row indices (first column),
%              column indices (second column) and values (third column).
%
% Note: this is a very general function to which direct calls are not
%       usually required - please use the (much friendlier) operator()
%       function.
%
% Note: the superoperator is returned in XYZ sparse format, which is 
%       different from Matlab's CSC format.
%
% ilya.kuprov@weizmann.ac.il
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=superop.m>

function A=superop(spin_system,opspec,side)

% Issue a recursive call if appropriate
if strcmp(side,'comm')
    A=superop(spin_system,opspec,'leftofcomm');
    B=superop(spin_system,opspec,'rightofcomm');
    B(:,3)=-B(:,3); A=[A; B]; return;
elseif strcmp(side,'acomm')
    A=[superop(spin_system,opspec,'left');
       superop(spin_system,opspec,'right')]; return;
end

% Validate the input
grumble(spin_system,opspec);

% Determine the relevant spins
active_spins=find(opspec);

% For unit operator use a shortcut
if isempty(active_spins)
    A=unit_oper(spin_system); 
    [rows,cols,vals]=find(A);
    A=[rows cols vals];
    return;
end

% Preallocate source state index
source=cell(1,numel(active_spins));

% Preallocate destination state index
destin=cell(1,numel(active_spins));

% Preallocate structure coefficients table
struct=cell(1,numel(active_spins));

% Loop over the relevant spins
for n=1:length(active_spins)
    
    % Get right and left product action tables for the current spin
    [pt_left,pt_right]=ist_product_table(spin_system.comp.mults(active_spins(n)));
    
    % Extract pages corresponding to the current state
    switch side
        case {'left','leftofcomm'}
            pt=squeeze(pt_left(opspec(active_spins(n))+1,:,:));
        case {'right','rightofcomm'}
            pt=squeeze(pt_right(opspec(active_spins(n))+1,:,:));
        otherwise
            error('invalid side specification.');
    end
    
    % Convert product action table to indices
    [destin{n},source{n},struct{n}]=find(pt);
    
    % Switch to 0 index for unit state
    source{n}=source{n}-1; destin{n}=destin{n}-1;
    
end
    
% Compute the structure coefficients for the relevant sub-algebra
from=source{1}; to=destin{1}; coeff=struct{1};
for n=2:numel(active_spins)
    from=[kron(from,ones(size(source{n},1),1)) ...
          kron(ones(size(from,1),1),source{n})];
    to=[kron(to,ones(size(destin{n},1),1)) ...
        kron(ones(size(to,1),1),destin{n})];
    coeff=kron(coeff,struct{n});
end

% Lift the basis columns corresponding to the relevant spins
basis_cols=spin_system.bas.basis(:,active_spins);

% For commutation superoperators remove commuting paths
if ismember(side,{'leftofcomm','rightofcomm'})
    kill_mask=(sum(from,2)==0)|(sum(to,2)==0);
    from(kill_mask,:)=[]; to(kill_mask,:)=[]; coeff(kill_mask,:)=[];
end

% Get the index array going
A=cell(size(from,1),1);

% Loop over source states
for n=1:size(from,1)
    
    % Retrieve the source subspace
    source_subsp_idx=true(size(basis_cols,1),1);
    for m=1:size(from,2)
        source_subsp_idx=and(source_subsp_idx,(basis_cols(:,m)==from(n,m)));
    end
    source_subsp=spin_system.bas.basis(source_subsp_idx,:);
    source_subsp_idx=find(source_subsp_idx);
    source_subsp(:,active_spins)=[];
    
    % Get source subspace dimension
    subsp_dim=size(source_subsp,1);
    
    % Skip empty subspaces
    if subsp_dim>0
    
        % Retrieve the destination subspace
        destin_subsp_idx=true(size(basis_cols,1),1);
        for m=1:size(to,2)
            destin_subsp_idx=and(destin_subsp_idx,(basis_cols(:,m)==to(n,m)));
        end
        destin_subsp=spin_system.bas.basis(destin_subsp_idx,:);
        destin_subsp_idx=find(destin_subsp_idx);
        destin_subsp(:,active_spins)=[];
        
        % Fill the operator
        if isequal(source_subsp,destin_subsp)
            
            % If the subspaces fully match, use the raw indices
            A{n}=[source_subsp_idx destin_subsp_idx coeff(n)*ones(subsp_dim,1)];
            
        else
        
            % Otherwise, use brute-force state-by-state matching
            [does_it_go_anywhere,where_it_goes_if_it_does]=ismember(source_subsp,destin_subsp,'rows');
            A{n}=[source_subsp_idx(does_it_go_anywhere)                           ...
                  destin_subsp_idx(where_it_goes_if_it_does(does_it_go_anywhere)) ...
                  coeff(n)*ones(nnz(does_it_go_anywhere),1)];
        
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

end

% Consistency enforcement
function grumble(spin_system,opspec)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function only supports sphten-liouv formalism.');
end
if (~isnumeric(opspec))||(~isrow(opspec))||any(mod(opspec,1)~=0)
    error('opspec must be a row vector of integers.');
end
if numel(opspec)~=spin_system.comp.nspins
    error('the number of elements in the opspec array must be equal to the number of spins.');
end
if any((opspec+1)>spin_system.comp.mults.^2)
    error('physically impossible state requested in opspec.');
end
end

% My philosophy, in essence, is the concept of man as a heroic being,
% with his own happiness as the moral purpose of his life, with pro-
% ductive achievement as his noblest activity, and reason as his only
% absolute.
%
% Ayn Rand

