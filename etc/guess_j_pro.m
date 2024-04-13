% Assigns J-couplings from literature values and Karplus curves. Syntax:
%
%            jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)
%
% Parameters:
%
%       aa_num   - a vector of amino acid numbers
%
%       aa_typ   - a cell array of amino acid types
%
%       pdb_id   - a cell array of PDB atom identifiers
%
%       coords   - a cell array of coordinate vectors
%
% Database nomenclature is as follows:
%
% 1. Atoms in the subgraph descriptor are listed alphabetically to make
%    the descriptors unique.
%
% 2. The four numbers refer to the bonding order, e.g. 1, 3, 2, 4 means
%    that the first atom in the descriptor is bonded to the third, which
%    is bonded to the second, which is bonded to the fourth. The coupling
%    in this case is between atom 1 and atom 4 in the descriptor.
%
% Note: these J-couplings should be considered approximate. For accurate
%       protein work you must supply your own J-couplings.
%
% Note: this is an auxiliary function that is called by protein.m protein
%       import module. Direct calls are discouraged.
%
% i.kuprov@soton.ac.uk
% ztw1e12@soton.ac.uk (Zenawi Welderufael)
% Andras_Boeszoermenyi@hms.harvard.edu
%
% <https://spindynamics.org/wiki/index.php?title=Guess_j_pro.m>

function jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)

% Check consistency
grumble(aa_num,aa_typ,pdb_id,coords);

% Preallocate the answer
jmatrix=cell(numel(coords),numel(coords));

% Number the atoms
numbers=1:numel(coords);

% Determine the connectivity graph
proxmatrix=false(numel(coords),numel(coords));
for n=1:numel(coords)
    for k=1:numel(coords)
        if (norm(coords{n}-coords{k},2)<1.60), proxmatrix(n,k)=1; end
    end
end

% Get all connected subgraphs up to size 2
subgraphs=dfpt(sparse(proxmatrix),2);

% Set generic coupling values
J_NH=-90.0; J_CN=-15.0; J_CH=140.0; J_CC=35.0;

% Spec all ordered connected pairs
pairs_database={...

% Backbone has specific numbers
'CA_CB'     , +34.9;   'CA_HA'     , +143.5;   'CA_N'      , -10.7;
'C_CA'      , +52.5;   'C_N'       ,  -14.4;   'H_N'       , -93.3;

% The rest is generic
'CA_HA2'    , J_CH;   'CA_HA3'    , J_CH;   'CB_CG'     , J_CC;   'CB_CG1'    , J_CC;
'CB_CG2'    , J_CC;   'CB_HB'     , J_CH;   'CB_HB1'    , J_CH;   'CB_HB2'    , J_CH;
'CB_HB3'    , J_CH;   'CD1_CE1'   , J_CC;   'CD1_CG'    , J_CC;   'CD1_CG1'   , J_CC;
'CD1_HD1'   , J_CH;   'CD1_HD11'  , J_CH;   'CD1_HD12'  , J_CH;   'CD1_HD13'  , J_CH;    'CD1_NE1'   , J_CN;
'CD2_CE2'   , J_CC;   'CD2_CG'    , J_CC;   'CD2_HD2'   , J_CH;   'CD2_HD21'  , J_CH;    'CD2_CE3'   , J_CC;
'CD2_HD22'  , J_CH;   'CD2_HD23'  , J_CH;   'CD2_NE2'   , J_CN;   'CD_CE'     , J_CC;    'CE2_NE1'   , J_CN;
'CD_CG'     , J_CC;   'CD_HD2'    , J_CH;   'CD_HD3'    , J_CH;   'CD_N'      , J_CN;    'HE1_NE1'   , J_NH;
'CD_NE2'    , J_CN;   'CE1_CZ'    , J_CC;   'CE1_HE1'   , J_CH;   'CE1_ND1'   , J_CN;    'CE2_CZ2'   , J_CC;
'CE1_NE2'   , J_CN;   'CE2_CZ'    , J_CC;   'CE2_HE2'   , J_CH;   'CE_HE1'    , J_CH;    'CE3_CZ3'   , J_CC;
'CE_HE2'    , J_CH;   'CE_HE3'    , J_CH;   'CG1_HG11'  , J_CH;   'CG1_HG12'  , J_CH;    'CE3_HE3'   , J_CH;
'CG1_HG13'  , J_CH;   'CG2_HG21'  , J_CH;   'CG2_HG22'  , J_CH;   'CG2_HG23'  , J_CH;    'CH2_CZ2'   , J_CC;
'CG_HG'     , J_CH;   'CH2_HH2'   , J_CH;   'CG_HG2'    , J_CH;   'CG_HG3'    , J_CH;    'CG_ND1'    , J_CN;
'CG_ND2'    , J_CN;   'CZ_HZ'     , J_CH;   'HD21_ND2'  , J_NH;   'HD22_ND2'  , J_NH;    'CZ3_HZ3'   , J_CH;
'HE21_NE2'  , J_NH;   'HE22_NE2'  , J_NH;   'HH22_NH2'  , J_NH;   'HH21_NH2'  , J_NH;    'CH2_CZ3'   , J_CC;
'HH12_NH1'  , J_NH;   'HH11_NH1'  , J_NH;   'CZ_NH2'    , J_CN;   'CZ_NH1'    , J_CN;    'CZ2_HZ2'   , J_CH;
'CZ_NE'     , J_CN;   'HE_NE'     , J_NH;   'CD_NE'     , J_CN;   'HZ3_NZ'    , J_NH;
'HZ2_NZ'    , J_NH;   'HZ1_NZ'    , J_NH;   'CE_NZ'     , J_CN;   'HD1_ND1'   , J_NH;
'H_H'       , 0.00;   'H1_N'      , J_NH;   'H2_N'      , J_NH;   'H3_N'      , J_NH};

% Loop over subgraphs
disp(' '); disp('############# ONE-BOND J-COUPLING SUMMARY ###############');
for n=1:size(subgraphs,1)
    
    % Extract labels, residues and numbers
    spin_labels=pdb_id(subgraphs(n,:));
    spin_numbers=numbers(subgraphs(n,:));
    spin_resnames=aa_typ(subgraphs(n,:));
    spin_resnums=aa_num(subgraphs(n,:));
    
    % Index the spins
    [spin_labels,index]=sort(spin_labels);
    spin_numbers=spin_numbers(index);
    spin_resnames=spin_resnames(index);
    spin_resnums=spin_resnums(index);
    
    % Consult the database
    if isscalar(spin_labels)
        
        % Inform the user about solitary spins
        disp(['WARNING: solitary spin detected, ' spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1}]);
        
    else
        
        % Build specification string
        spec_string=pairs_database(strcmp([spin_labels{1} '_' spin_labels{2}],pairs_database(:,1)),:);
        
        % Assign the coupling
        if isempty(spec_string)
            
            % Complain if nothing found
            disp(['WARNING: unknown atom pair or unphysical proximity, ' spin_labels{1} '_' spin_labels{2} ' in residues ' spin_resnames{1}...
                  '(' num2str(spin_resnums(1)) '), ' spin_resnames{2} '(' num2str(spin_resnums(2)) ')']);
        
        elseif size(spec_string,1)>1
            
            % Complain if multiple records found
            error([spec_string{1} ' coupling has been specified multiple times in the coupling database.']);
            
        else
    
            % Write the array
            jmatrix{spin_numbers(1),spin_numbers(2)}=spec_string{2};
    
            % Report to the user
            disp(['Estimated one-bond J-coupling from ' pad([spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1}],19) ' to ' ...
                                                        pad([spin_resnames{2} '(' num2str(spin_resnums(2)) '):' spin_labels{2}],19) ' ' ...
                                                        pad(num2str(spec_string{2},'%5.1f'),6,'left') ' Hz']);
            
        end
        
    end
                                            
end

% Get all connected subgraphs up to size 3
subgraphs=dfpt(sparse(proxmatrix),3);

% Set generic coupling values
J_CCC=-1.2;  % Zenawi's DFT
J_CCH=0;     % Literature reference
J_CCN=7.0;   % Zenawi's DFT
J_HCH=-12.0; % Literature reference
J_HNH=-1.5;  % Zenawi's DFT
J_NCN=0;     % Literature reference
J_NCH=0;     % Literature reference
J_CNH=0;     % Literature reference
J_CNC=0;     % Literature reference

% Spec all ordered connected triples (sorted atoms, bonding order, coupling from atom 1 to atom 3)
triples_database={...
    
% Backbone has specific numbers
'C_CA_CB',       1, 2, 3,  -0.4;    'C_CA_HA',       1, 2, 3, -4.5;    'C_CA_HA2',      1, 2, 3, -4.5;    'C_CA_HA3',      1, 2, 3,  -4.5;
'C_CA_N',        2, 1, 3,  -7.7;    'C_CD_N',        1, 3, 2,  0.0;    'C_H_N',         1, 3, 2,  4.0;    'CA_CB_HA',      2, 1, 3,  -4.6;
'CA_CB_N',       2, 1, 3,  -0.5;    'CA_H_N',        1, 3, 2,  2.0;    'CA_HA_N',       2, 1, 3,  0.5;    'CA_HA2_HA3',    2, 1, 3, -12.0;
'CA_HA2_N',      2, 1, 3,   0.5;    'CA_HA3_N',      2, 1, 3,  0.5;    'CA_CD_N',       1, 3, 2,  0.0;

% Everything else is generic
'CA_CB_CG',      1, 2, 3, -1.2;    'CA_CB_CG1',     1, 2, 3, -1.2;    'CA_CB_CG2',     1, 2, 3, J_CCC;   'CA_CB_HB',      1, 2, 3, J_CCH;
'CA_CB_HB1',     1, 2, 3, -7.5;    'CA_CB_HB2',     1, 2, 3, -4.2;    'CA_CB_HB3',     1, 2, 3, -2.2;    'CB_CD_CG',      1, 3, 2, J_CCC;
'CB_CD1_CG',     1, 3, 2, J_CCC;   'CB_CD1_CG1',    1, 3, 2, J_CCC;   'CB_CD2_CG',     1, 3, 2, J_CCC;   'CB_CG_HB2',     2, 1, 3, J_CCH;
'CB_CG_HB3',     2, 1, 3, J_CCH;   'CB_CG_HG',      1, 2, 3, J_CCH;   'CB_CG_HG2',     1, 2, 3, J_CCH;   'CB_CG_HG3',     1, 2, 3, J_CCH;
'CB_CG_ND1',     1, 2, 3, J_CCN;   'CB_CG_ND2',     1, 2, 3, J_CCN;   'CB_CG1_CG2',    2, 1, 3, J_CCC;   'CB_CG1_HB',     2, 1, 3, J_CCH;
'CB_CG1_HG11',   1, 2, 3, J_CCH;   'CB_CG1_HG12',   1, 2, 3, J_CCH;   'CB_CG1_HG13',   1, 2, 3, J_CCH;   'CB_CG2_HB',     2, 1, 3, J_CCH;
'CB_CG2_HG21',   1, 2, 3, J_CCH;   'CB_CG2_HG22',   1, 2, 3, J_CCH;   'CB_CG2_HG23',   1, 2, 3, J_CCH;   'CB_HB1_HB2',    2, 1, 3, J_HCH;
'CB_HB1_HB3',    2, 1, 3, J_HCH;   'CB_HB2_HB3',    2, 1, 3, J_HCH;   'CD_CE_CG',      2, 1, 3, J_CCC;   'CD_CE_HD2',     2, 1, 3, J_CCH;
'CD_CE_HD3',     2, 1, 3, J_CCH;   'CD_CE_HE2',     1, 2, 3, J_CCH;   'CD_CE_HE3',     1, 2, 3, J_CCH;   'CD_CG_HD2',     2, 1, 3, J_CCH;
'CD_CG_HD3',     2, 1, 3, J_CCH;   'CD_CG_HG2',     1, 2, 3, J_CCH;   'CD_CG_HG3',     1, 2, 3, J_CCH;   'CD_CG_N',       2, 1, 3, J_CCN;
'CD_CG_NE2',     2, 1, 3, J_CCN;   'CD_HD2_HD3',    2, 1, 3, J_HCH;   'CD_HD2_N',      2, 1, 3, J_NCH;   'CD_HD3_N',      2, 1, 3, J_NCH;
'CD_HE21_NE2',   1, 3, 2, J_CNH;   'CD_HE22_NE2',   1, 3, 2, J_CNH;   'CD1_CD2_CG',    1, 3, 2, J_CCC;   'CD1_CE1_CG',    2, 1, 3, J_CCC;
'CD1_CE1_CZ',    1, 2, 3, J_CCC;   'CD1_CE1_HD1',   2, 1, 3, J_CCH;   'CD1_CE1_HE1',   1, 2, 3, J_CCH;   
'CD1_CG_HD11',   2, 1, 3, J_CCH;   'CD1_CG_HD12',   2, 1, 3, J_CCH;   'CD1_CG_HD13',   2, 1, 3, J_CCH;   'CD1_CG_HG',     1, 2, 3, J_CCH;
'CD1_CG1_HD11',  2, 1, 3, J_CCH;   'CD1_CG1_HD12',  2, 1, 3, J_CCH;   'CD1_CG1_HD13',  2, 1, 3, J_CCH;   'CD1_CG1_HG12',  1, 2, 3, J_CCH;
'CD1_CG1_HG13',  1, 2, 3, J_CCH;   'CD1_HD11_HD12', 2, 1, 3, J_HCH;   'CD1_HD11_HD13', 2, 1, 3, J_HCH;   'CD1_HD12_HD13', 2, 1, 3, J_HCH;    
'CD2_CE1_NE2',   1, 3, 2, J_CNC;   'CD2_CE2_CG',    2, 1, 3, J_CCC;   'CD2_CE2_CZ',    1, 2, 3, J_CCC;   'CD2_CE2_HD2',   2, 1, 3, J_CCH;    
'CD2_CE2_HE2',   1, 2, 3, J_CCH;   'CD2_CG_HD2',    2, 1, 3, J_CCH;   'CD2_CG_HD21',   2, 1, 3, J_CCH;   'CD2_CG_HD22',   2, 1, 3, J_CCH;    
'CD2_CG_HD23',   2, 1, 3, J_CCH;   'CD2_CG_HG',     1, 2, 3, J_CCH;   'CD2_CG_ND1',    1, 2, 3, J_CCN;   'CD2_CG_NE2',    2, 1, 3, J_CCN;    
'CD2_HD2_NE2',   2, 1, 3, J_NCH;   'CD2_HD21_HD22', 2, 1, 3, J_HCH;   'CD2_HD21_HD23', 2, 1, 3, J_HCH;   'CD2_HD22_HD23', 2, 1, 3, J_HCH;    
'CE_HE1_HE2',    2, 1, 3, J_HCH;   'CE_HE1_HE3',    2, 1, 3, J_HCH;   'CE_HE2_HE3',    2, 1, 3, J_HCH;   'CE1_CE2_CZ',    1, 3, 2, J_CCC;    
'CE1_CG_ND1',    1, 3, 2, J_CNC;   'CE1_CZ_HE1',    2, 1, 3, J_CCH;   'CE1_CZ_HZ',     1, 2, 3, J_CCH;   'CE1_HE1_ND1',   2, 1, 3, J_NCH;    
'CE1_HE1_NE2',   2, 1, 3, J_NCH;   'CE1_ND1_NE2',   2, 1, 3, J_NCN;   'CE2_CZ_HE2',    2, 1, 3, J_CCH;   'CE2_CZ_HZ',     1, 2, 3, J_CCH;   
'CG_HD21_ND2',   1, 3, 2, J_CNH;   'CG_HD22_ND2',   1, 3, 2, J_CNH;   'CG_HG2_HG3',    2, 1, 3, J_HCH;   'CG1_HG11_HG12', 2, 1, 3, J_HCH;    
'CG1_HG11_HG13', 2, 1, 3, J_HCH;   'CG1_HG12_HG13', 2, 1, 3, J_HCH;   'CG2_HG21_HG22', 2, 1, 3, J_HCH;   'CG2_HG21_HG23', 2, 1, 3, J_HCH;    
'CG2_HG22_HG23', 2, 1, 3, J_HCH;   'HD21_HD22_ND2', 1, 3, 2, J_HNH;   'HE21_HE22_NE2', 1, 3, 2, J_HNH;   'HH21_HH22_NH2', 1, 3, 2, J_HNH;    
'HH11_HH12_NH1', 1, 3, 2, J_HNH;   'CZ_HH22_NH2',   1, 3, 2, J_CNH;   'CZ_HH21_NH2',   1, 3, 2, J_CNH;   'CZ_NH1_NH2',    2, 1, 3, J_NCN;    
'CZ_HH12_NH1',   1, 3, 2, J_CNH;   'CZ_HH11_NH1',   1, 3, 2, J_CNH;   'CZ_NE_NH2',     2, 1, 3, J_NCN;   'CZ_NE_NH1',     2, 1, 3, J_NCN;     
'CZ_HE_NE',      1, 3, 2, J_CNH;   'CD_CZ_NE',      1, 3, 2, J_CNC;   'CD_HE_NE',      1, 3, 2, J_CNH;   'CD_HD3_NE',     2, 1, 3, J_NCH;     
'CD_HD2_NE',     2, 1, 3, J_NCH;   'CD_CG_NE',      2, 1, 3, J_CCN;   'HZ2_HZ3_NZ'     1, 3, 2, J_HNH;   'HZ1_HZ3_NZ',    1, 3, 2, J_HNH;    
'HZ1_HZ2_NZ',    1, 3, 2, J_HNH;   'CE_HZ3_NZ',     1, 3, 2, J_CNH;   'CE_HZ2_NZ',     1, 3, 2, J_CNH;   'CE_HZ1_NZ',     1, 3, 2, J_CNH;    
'CE_HE3_NZ',     2, 1, 3, J_NCH;   'CE_HE2_NZ',     2, 1, 3, J_NCH;   'CD_CE_NZ',      1, 2, 3, J_CCN;   'CE1_HD1_ND1',   1, 3, 2, J_CNH;    
'CG_HD1_ND1',    1, 3, 2, J_CNH;   'H_H_N',         1, 3, 2, J_HNH;   'H2_H3_N',       1, 3, 2, J_HNH;   'H1_H2_N',       1, 3, 2, J_HNH;
'H1_H3_N',       1, 3, 2, J_HNH;   'CA_H1_N',       1, 3, 2, J_CNH;   'CA_H2_N',       1, 3, 2, J_CNH;   'CA_H3_N',       1, 3, 2, J_CNH;

% Tryptophan side chain (GIAO DFT M06/cc-pVTZ data in PCM water)
'CH2_CZ3_HH2',   2, 1, 3, -1.82;   'CH2_CZ3_HZ3',   1, 2, 3, -1.13;   'CH2_CZ2_HH2',   2, 1, 3,  -2.82;   'CH2_CZ2_HZ2',   1, 2, 3, -3.81;
'CH2_CZ2_CZ3',   2, 1, 3, -3.59;   'CE3_CZ3_HZ3',   1, 2, 3, -2.22;   'CE3_CZ3_HE3',   2, 1, 3,  -3.50;   'CE3_CH2_CZ3',   1, 3, 2, -3.38;
'CE2_CZ2_HZ2',   1, 2, 3, -3.47;   'CE2_CH2_CZ2',   1, 3, 2, -1.27;   'CE2_HE1_NE1',   1, 3, 2,   2.36;   'CE2_CZ2_NE1',   2, 1, 3,  0.72;
'CD2_CE3_HE3',   1, 2, 3, -1.35;   'CD2_CE3_CZ3',   1, 2, 3, -2.69;   'CD2_CE2_CZ2',   1, 2, 3,   0.26;   'CD2_CE2_CE3',   2, 1, 3,  0.61;
'CD2_CE2_NE1',   1, 2, 3,  3.55;   'CD1_HE1_NE1',   1, 3, 2,  3.38;   'CD1_HD1_NE1',   2, 1, 3,   2.72;   'CD1_CE2_NE1',   1, 3, 2,  5.55;
'CD2_CE3_CG',    3, 1, 2,  1.08;   'CD1_CG_NE1',    2, 1, 3,  2.30};

% Collision detection exceptions (HIS and PRO, where two-bond and three-bond couplings collide)
allowed_collision_list={'H_H_N'};

% Loop over subgraphs
disp(' '); disp('############# TWO-BOND J-COUPLING SUMMARY ###############');
for n=1:size(subgraphs,1)
    
    % Extract labels, residues and numbers
    spin_labels=pdb_id(subgraphs(n,:));
    spin_numbers=numbers(subgraphs(n,:));
    spin_resnames=aa_typ(subgraphs(n,:));
    spin_resnums=aa_num(subgraphs(n,:));
    
    % Index the spins
    [spin_labels,index]=sort(spin_labels);
    spin_numbers=spin_numbers(index);
    spin_resnames=spin_resnames(index);
    spin_resnums=spin_resnums(index);
    
    % Consult the database
    if isscalar(spin_labels)
        
        % Inform the user about solitary spins
        disp(['WARNING: solitary spin detected, ' spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1}]);
        
    elseif numel(spin_labels)==2
        
        % Inform the user about isolated pairs
        disp(['WARNING: isolated spin pair detected, ' spin_resnames{1} ':' spin_labels{1} ' and ' spin_resnames{2} ':' spin_labels{2}]);
        
    else
        
        % Build the specification string
        spec_string=triples_database(strcmp([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3}],triples_database(:,1)),:);
    
        % Assign the coupling
        if isempty(spec_string)
            
            % Inform the user about missing data
            disp(['WARNING: unknown atom triple or unphysical proximity, ' spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} ' in amino acids '...
                                                                           spin_resnames{1} '(' num2str(spin_resnums(1)) '), '...
                                                                           spin_resnames{2} '(' num2str(spin_resnums(2)) '), '...
                                                                           spin_resnames{3} '(' num2str(spin_resnums(3)) ')']);
        
        elseif size(spec_string,1)>1
            
            % Complain if multiple records found
            error([spec_string{1} ' coupling has been specified multiple times in the coupling database.']);    
        
        else
            
            % Extract spin numbers
            spin_a=spin_numbers(spec_string{2}); spin_b=spin_numbers(spec_string{3}); spin_c=spin_numbers(spec_string{4});
            
            % Match bonding pattern
            if proxmatrix(spin_a,spin_b)&&proxmatrix(spin_b,spin_c)
                
                % Check for collisions
                if (~isempty(jmatrix{spin_a,spin_c}))&&(~ismember([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3}],allowed_collision_list))
                    
                    % Complain and bomb out
                    error(['Atom triple: ' spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} ' in amino acids '...
                                           spin_resnames{1} '(' num2str(spin_resnums(1)) '), '...
                                           spin_resnames{2} '(' num2str(spin_resnums(2)) '), '...
                                           spin_resnames{3} '(' num2str(spin_resnums(3)) ' collides with other data.']);
                    
                else
                    
                    % Write the array
                    jmatrix{spin_a,spin_c}=spec_string{5};
                    
                end
                
                % Inform the user
                disp(['Estimated two-bond J-coupling from ' pad([aa_typ{spin_a} '(' num2str(aa_num(spin_a)) '):' pdb_id{spin_a}],19) ' to '...
                                                            pad([aa_typ{spin_c} '(' num2str(aa_num(spin_c)) '):' pdb_id{spin_c}],19) ' '...
                                                            pad(num2str(jmatrix{spin_a,spin_c},'%5.1f'),6,'left') ' Hz']);
                
            end
            
        end
        
    end

end

% Get all connected subgraphs up to size 4
subgraphs=dfpt(sparse(proxmatrix),4);

% Remove T-shaped subgraphs
kill_mask=false(size(subgraphs,1),1);
for n=1:size(subgraphs,1)
    if any(sum(proxmatrix(subgraphs(n,:),subgraphs(n,:)))==4)
        kill_mask(n)=true();
    end
end
subgraphs(kill_mask,:)=[];

% Generic Karplus coefficients
J_CCCC_AL=[4.48  0.18 -0.57]; % Aliphatic C-C-C-C from http://dx.doi.org/10.1016/j.carres.2007.02.023
J_CCCH_AL=[3.83 -0.90  3.81]; % Aliphatic C-C-C-H from http://onlinelibrary.wiley.com/doi/10.1002/anie.198204491/pdf
J_HCCH_AL=[4.22 -0.50  4.50]; % Aliphatic H-C-C-H from http://pubs.acs.org/doi/abs/10.1021/ja00901a059
J_CCCC_AR=[0.00  0.00  5.00]; % Aromatic  C-C-C-C from [literature reference]
J_CCCH_AR=[0.00  0.00  8.20]; % Aromatic  C-C-C-H from [literature reference]
J_HCCH_AR=[0.00  0.00  10.2]; % Aromatic  H-C-C-H from [literature reference]

% List of ordered connected quads (sorted atoms, bonding order, Karplus coefficients)
quads_database={...

% Backbone and its vicinity
'CA_CB_CG1_HA',      4, 1, 2, 3, [7.10 -1.00 0.70];    % Vuister, G. W. "J-couplings. Measurement and Usage in Structure Determination." [chi1 side-chain]
'CA_CB_CG2_HA',      4, 1, 2, 3, [7.10 -1.00 0.70];    % Vuister, G. W. "J-couplings. Measurement and Usage in Structure Determination." [chi1 side-chain]
'CA_CB_CG_H',        1, 2, 3, 4, [4.03 -0.87 4.50];    % http://dx.doi.org/10.1002/mrc.1260280513    chi2 ??? need to look at it later. 
'CA_CB_CG1_HA',      4, 1, 2, 3, [10.2 -1.30 0.20];    % http://dx.doi.org/10.1021/bi00585a003              (checked)
'CA_CB_CG2_HA',      4, 1, 2, 3, [10.2 -1.30 0.20];    % http://dx.doi.org/10.1021/bi00585a003              (checked)
'CA_CB_CG1_HA',      4, 1, 2, 3, [10.2 -1.30 0.20];    % A. DeMarco and M. Llinas, Biochemistry 18, 3846 (1979)    [chi1 side-chain]       (checked)
'CA_CB_CG2_HA',      4, 1, 2, 3, [10.2 -1.30 0.20];    % A. DeMarco and M. Llinas, Biochemistry 18, 3846 (1979)    [chi1 side-chain]       (checked)
'CA_CB_CG_HA',       4, 1, 2, 3, [10.2 -1.30 0.20];    % A. DeMarco and M. Llinas, Biochemistry 18, 3846 (1979)    [chi1 side-chain]       (checked)
'C_CA_CB_CG',        1, 2, 3, 4, [3.42 -0.59 0.17];    % http://pubs.acs.org/doi/pdf/10.1021/ja029972s             [chi1 side-chain]       (checked)
'C_CA_CB_CG1',       1, 2, 3, 4, [3.42 -0.59 0.17];    % http://pubs.acs.org/doi/pdf/10.1021/ja029972s             [chi1 side-chain]       (checked)
'C_CA_CB_CG2',       1, 2, 3, 4, [3.30 -0.51 0.04];    % http://pubs.acs.org/doi/pdf/10.1021/ja029972s             [chi1 side-chain]       (checked)
'CA_CB_CG_N',        4, 1, 2, 3, [2.64 0.26 -0.22];    % http://pubs.acs.org/doi/pdf/10.1021/ja029972s             [chi1 side-chain]       (checked)
'CA_CB_CG1_N',       4, 1, 2, 3, [2.22 0.25 -0.06];    % http://pubs.acs.org/doi/pdf/10.1021/ja029972s             [chi1 side-chain]       (checked) 
'CA_CB_CG2_N',       4, 1, 2, 3, [2.01 0.21 -0.12];    % http://pubs.acs.org/doi/pdf/10.1021/ja029972s             [chi1 side-chain]       (checked)
'C_CA_CB_HB3',       1, 2, 3, 4, [7.2  -2.04  0.60];   % http://pubs.acs.org/doi/pdf/10.1021/ja00528a004           [chi1 side-chain]       (checked)
'C_CA_CB_HB2',       1, 2, 3, 4, [7.2  -2.04  0.60];   % http://pubs.acs.org/doi/pdf/10.1021/ja00528a004           [chi1 side-chain]       (checked)
'C_CA_CB_HB1',       1, 2, 3, 4, [7.2  -2.04  0.60];   % http://pubs.acs.org/doi/pdf/10.1021/ja00528a004           [chi1 side-chain]       (checked)
'C_CA_CB_HB',        1, 2, 3, 4, [7.2  -2.04  0.60];   % http://pubs.acs.org/doi/pdf/10.1021/ja00528a004           [chi1 side-chain]       (checked)
'CA_CB_HA_HB',       3, 1, 2, 4, [ 9.5 -1.60  1.80];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)
'CA_CB_HA_HB1',      3, 1, 2, 4, [ 9.5 -1.60  1.80];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)
'CA_CB_HA_HB2',      3, 1, 2, 4, [ 9.5 -1.60  1.80];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)  
'CA_CB_HA_HB3',      3, 1, 2, 4, [ 9.5 -1.60  1.80];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)
'CA_CB_HB1_N',       4, 1, 2, 3, [-4.40  1.20 0.10];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)
'CA_CB_HB2_N',       4, 1, 2, 3, [-4.40  1.20 0.10];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)
'CA_CB_HB3_N',       4, 1, 2, 3, [-4.40  1.20 0.10];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)
'CA_CB_HB_N',        4, 1, 2, 3, [-4.40  1.20 0.10];   % http://dx.doi.org/10.1006/jmbi.1998.1911                  [chi1 side-chain]       (checked)
'CA_H_HA3_N',        2, 4, 1, 3, [6.51 -1.76 1.60];    % http://dx.doi.org/10.1021/ja00070a024                     [backbone phi]          (checked)
'CA_H_HA_N',         2, 4, 1, 3, [6.51 -1.76 1.60];    % http://dx.doi.org/10.1021/ja00070a024                     [backbone phi]          (checked)
'CA_H_HA2_N',        2, 4, 1, 3, [6.51 -1.76 1.60];    % http://dx.doi.org/10.1021/ja00070a024                     [backbone phi]          (checked)
'CA_CB_H_N',         3, 4, 1, 2, [3.06 -0.074 0.13];   % http://pubs.acs.org/doi/pdf/10.1021/ja970067v             [backbone phi]          (checked)
'C_CA_H_N',          1, 2, 4, 3, [4.29 -1.01 0.00];    % http://pubs.acs.org/doi/pdf/10.1021/ja001798p             [backbone phi]          (checked)
'C_CA_HA2_N',        1, 4, 2, 3, [3.72 -2.18 1.28];    % http://pubs.acs.org/doi/pdf/10.1021/ja001798p             [backbone phi]          (checked)                                                         
'C_CA_HA3_N',        1, 4, 2, 3, [3.72 -2.18 1.28];    % http://pubs.acs.org/doi/pdf/10.1021/ja001798p             [backbone phi]          (checked)
'C_CA_HA_N',         1, 4, 2, 3, [3.72 -2.18 1.28];    % http://pubs.acs.org/doi/pdf/10.1021/ja001798p             [backbone phi]          (checked)
'C_CA_CB_N',         1, 4, 2, 3, [1.59 -0.67 0.27];    % http://spin.niddk.nih.gov/bax/lit/508/244.pdf             [backbone phi]          (checked)
'C_C_CA_N',          1, 4, 3, 2, [1.33 -0.88 0.06];    % http://pubs.acs.org/doi/pdf/10.1021/ja9616239             [backbone phi]          (checked)
'C_CB_CA_N',         1, 4, 3, 2, [1.59 -0.67 0.27];    % http://spin.niddk.nih.gov/bax/lit/508/244.pdf             [backbone phi]          (checked)
'C_CA_HA2_N',        4, 1, 2, 3, [-0.88 -0.61 -0.27];  % http://pubs.acs.org/doi/pdf/10.1021/ja00111a021           [backbone psi]          (checked)
'C_CA_HA3_N',        4, 1, 2, 3, [-0.88 -0.61 -0.27];  % http://pubs.acs.org/doi/pdf/10.1021/ja00111a021           [backbone psi]          (checked)
'C_CA_HA_N',         4, 1, 2, 3, [-0.88 -0.61 -0.27];  % http://pubs.acs.org/doi/pdf/10.1021/ja00111a021           [backbone psi]          (checked)
'C_CA_N_N',          3, 2, 1, 4, [0.0 0.0 0.0];        % To be filled in  
'C_CA_N_N',          4, 2, 1, 3, [0.0 0.0 0.0];        % To be filled in
'C_CB_CA_N',         4, 1, 3, 2, [0.0 0.0 0.0];        % To be filled in
'C_CA_CA_N',         2, 1, 4, 3, [0.0 0.0 0.0];        % To be filled in

% Histidine ring (GIAO DFT M06/cc-pVTZ)
'CD2_CG_HD2_ND1',    3, 1, 2, 4, [0.0 0.0 2.8];   'CD2_CG_ND1_NE2',    4, 1, 2, 3, [0.0 0.0 0.0];
'CB_CG_HB2_ND1',     3, 1, 2, 4, [0.0 0.0 0.0];   'CB_CG_HB3_ND1',     3, 1, 2, 4, [0.0 0.0 0.0];   
'CB_CE1_CG_ND1',     1, 3, 4, 2, [0.0 0.0 2.3];   'CD2_CE1_HE1_NE2',   3, 2, 4, 1, [0.0 0.0 4.9];
'CD2_CE1_HD2_NE2',   3, 1, 4, 2, [0.0 0.0 6.6];   'CD2_CE1_CG_NE2',    3, 1, 4, 2, [0.0 0.0 1.0];
'CD2_CE1_ND1_NE2',   1, 4, 2, 3, [0.0 0.0 0.0];   'CD2_CE1_CG_ND1',    1, 3, 4, 2, [0.0 0.0 3.4];
'CE1_CG_HE1_ND1',    3, 1, 4, 2, [0.0 0.0 7.0];   'CE1_CG_ND1_NE2',    2, 3, 1, 4, [0.0 0.0 1.4];
'CA_CB_CG_ND1',      1, 2, 3, 4, [0.0 0.0 1.0];   'CB_CD2_CG_NE2',     1, 3, 2, 4, [0.0 0.0 0.0];
'CE1_HD1_HE1_ND1',   2, 4, 1, 3, [0.0 0.0 3.9];   'CE1_HD1_ND1_NE2',   2, 3, 1, 4, [0.0 0.0 3.3];
'CD2_CG_HD1_ND1',    3, 4, 2, 1, [0.0 0.0 6.1];   'CB_CG_HD1_ND1',     3, 4, 2, 1, [0.0 0.0 1.0];

% Proline ring (GIAO DFT M06/cc-pVTZ)
'CA_CD_HA_N',        3, 1, 4, 2, [0.0 0.0 0.0];   'CA_CD_HD2_N',       1, 4, 2, 3, [0.0 0.0 2.5];
'CA_CD_HD3_N',       1, 4, 2, 3, [0.0 0.0 2.5];   'CA_CD_CG_N',        1, 4, 2, 3, [0.0 0.0 0.0];
'CD_CG_HG2_N',       3, 2, 1, 4, [0.0 0.0 0.0];   'CD_CG_HG3_N',       3, 2, 1, 4, [0.0 0.0 0.0];  
'C_CD_CG_N',         1, 4, 2, 3, [0.0 0.0 2.5];   'C_CA_CD_N',         1, 2, 4, 3, [0.0 0.0 1.1];
'C_CD_HD2_N',        1, 4, 2, 3, [0.0 0.0 1.0];   'C_CD_HD3_N',        1, 4, 2, 3, [0.0 0.0 1.0];
'CA_CB_CD_N',        2, 1, 4, 3, [0.0 0.0 0.0];   'CB_CD_CG_N',        4, 2, 3, 1, [0.0 0.0 1.1];  

% Generic aliphatics
'CA_CB_CD1_CG',      1, 2, 4, 3, J_CCCC_AL;       'CA_CB_CD1_CG1',     1, 2, 4, 3, J_CCCC_AL; 
'CA_CB_CD2_CG',      1, 2, 4, 3, J_CCCC_AL;       'CA_CB_CD_CG',       1, 2, 4, 3, J_CCCC_AL;
'CA_CB_CG2_HG21',    1, 2, 3, 4, J_CCCH_AL;       'CA_CB_CG_HG3',      1, 2, 3, 4, J_CCCH_AL;
'CA_CB_CG1_HG11',    1, 2, 3, 4, J_CCCH_AL;       'CA_CB_CG1_HG12',    1, 2, 3, 4, J_CCCH_AL;
'CA_CB_CG1_HG13',    1, 2, 3, 4, J_CCCH_AL;       'CD2_CG_HD21_HG',    4, 2, 1, 3, J_HCCH_AL;
'CA_CB_CG2_HG22',    1, 2, 3, 4, J_CCCH_AL;       'CA_CB_CG2_HG23',    1, 2, 3, 4, J_CCCH_AL;
'CA_CB_CG_HG',       1, 2, 3, 4, J_CCCH_AL;       'CA_CB_CG_HG2',      1, 2, 3, 4, J_CCCH_AL;
'CD_CG_HD3_HG2',     4, 2, 1, 3, J_HCCH_AL;       'CB_CD1_CG1_HD11',   1, 3, 2, 4, J_CCCH_AL;
'CB_CD1_CG1_CG2',    2, 3, 1, 4, J_CCCC_AL;       'CB_CD1_CG1_HB',     4, 1, 3, 2, J_CCCH_AL;     
'CB_CD1_CG1_HD12',   1, 3, 2, 4, J_CCCH_AL;       'CB_CD1_CG1_HD13',   1, 3, 2, 4, J_CCCH_AL;
'CB_CD1_CG_HB2',     4, 1, 3, 2, J_CCCH_AL;       'CB_CD1_CG_HB3',     4, 1, 3, 2, J_CCCH_AL;
'CB_CD1_CG_HD1',     1, 2, 3, 4, J_CCCH_AL;       'CB_CD1_CG_HD11',    1, 2, 3, 4, J_CCCH_AL;
'CB_CD1_CG_HD12',    1, 2, 3, 4, J_CCCH_AL;       'CB_CD1_CG_HD13',    1, 2, 3, 4, J_CCCH_AL;
'CB_CD2_CE2_CG',     1, 4, 2, 3, J_CCCC_AL;       'CB_CD2_CG_HB2',     4, 1, 3, 2, J_CCCH_AL;
'CB_CD2_CG_HB3',     4, 1, 3, 2, J_CCCH_AL;       'CB_CD2_CG_HD2',     1, 3, 2, 4, J_CCCH_AL;
'CB_CD2_CG_HD21',    1, 3, 2, 4, J_CCCH_AL;       'CB_CD2_CG_HD22',    1, 3, 2, 4, J_CCCH_AL;
'CB_CD2_CG_HD23',    1, 3, 2, 4, J_CCCH_AL;       'CB_CD_CG_HD3',      1, 3, 2, 4, J_CCCH_AL;
'CB_CD_CE_CG',       1, 4, 2, 3, J_CCCC_AL;       'CB_CD_CG_HB2',      4, 1, 3, 2, J_CCCH_AL;
'CB_CD_CG_HB3',      4, 1, 3, 2, J_CCCH_AL;       'CB_CD_CG_HD2',      1, 3, 2, 4, J_CCCH_AL;
'CB_CG1_CG2_HG11',   4, 2, 1, 3, J_CCCH_AL;       'CB_CG1_CG2_HG12',   4, 2, 1, 3, J_CCCH_AL;
'CB_CG1_CG2_HG13',   4, 2, 1, 3, J_CCCH_AL;       'CB_CG1_CG2_HG21',   4, 3, 1, 2, J_CCCH_AL;
'CB_CG1_CG2_HG22',   4, 3, 1, 2, J_CCCH_AL;       'CB_CG1_CG2_HG23',   4, 3, 1, 2, J_CCCH_AL;
'CB_CG1_HB_HG11',    4, 2, 1, 3, J_HCCH_AL;       'CB_CG1_HB_HG12',    4, 2, 1, 3, J_HCCH_AL;
'CB_CG1_HB_HG13',    4, 2, 1, 3, J_HCCH_AL;       'CB_CG2_HB_HG21',    4, 2, 1, 3, J_HCCH_AL;
'CB_CG2_HB_HG22',    4, 2, 1, 3, J_HCCH_AL;       'CB_CG2_HB_HG23',    4, 2, 1, 3, J_HCCH_AL;
'CB_CG_HB2_HG',      4, 2, 1, 3, J_HCCH_AL;       'CB_CG_HB2_HG2',     4, 2, 1, 3, J_HCCH_AL;
'CB_CG_HB2_HG3',     4, 2, 1, 3, J_HCCH_AL;       'CB_CG_HB3_HG',      4, 2, 1, 3, J_HCCH_AL;
'CB_CG_HB3_HG2',     4, 2, 1, 3, J_HCCH_AL;       'CB_CG_HB3_HG3',     4, 2, 1, 3, J_HCCH_AL;
'CD1_CD2_CG_HD1',    4, 1, 3, 2, J_CCCH_AL;       'CD1_CD2_CG_HD11',   4, 1, 3, 2, J_CCCH_AL;
'CD1_CD2_CG_HD12',   4, 1, 3, 2, J_CCCH_AL;       'CD1_CD2_CG_HD13',   4, 1, 3, 2, J_CCCH_AL;
'CD1_CD2_CG_HD2',    4, 2, 3, 1, J_CCCH_AL;       'CD1_CD2_CG_HD21',   4, 2, 3, 1, J_CCCH_AL;
'CD1_CD2_CG_HD22',   4, 2, 3, 1, J_CCCH_AL;       'CD1_CD2_CG_HD23',   4, 2, 3, 1, J_CCCH_AL;
'CD1_CD2_CE1_CG',    3, 1, 4, 2, J_CCCC_AL;       'CD1_CD2_CE2_CG',    1, 4, 2, 3, J_CCCC_AL;
'CD1_CG1_HD11_HG12', 4, 2, 1, 3, J_HCCH_AL;       'CD1_CG1_HD11_HG13', 4, 2, 1, 3, J_HCCH_AL;
'CD1_CG1_HD12_HG12', 4, 2, 1, 3, J_HCCH_AL;       'CD1_CG1_HD12_HG13', 4, 2, 1, 3, J_HCCH_AL;
'CD1_CG1_HD13_HG12', 4, 2, 1, 3, J_HCCH_AL;       'CD1_CG1_HD13_HG13', 4, 2, 1, 3, J_HCCH_AL;
'CD1_CG_HD11_HG',    4, 2, 1, 3, J_HCCH_AL;       'CD1_CG_HD12_HG',    4, 2, 1, 3, J_HCCH_AL;
'CD1_CG_HD13_HG',    4, 2, 1, 3, J_HCCH_AL;       'CD_CG_HD3_HG3',     4, 2, 1, 3, J_HCCH_AL;
'CD2_CG_HD22_HG',    4, 2, 1, 3, J_HCCH_AL;       'CD2_CG_HD23_HG',    4, 2, 1, 3, J_HCCH_AL;
'CD_CE_CG_HE2',      4, 2, 1, 3, J_CCCH_AL;       'CD_CE_CG_HE3',      4, 2, 1, 3, J_CCCH_AL;
'CD_CE_CG_HG2',      4, 3, 1, 2, J_CCCH_AL;       'CD_CE_CG_HG3',      4, 3, 1, 2, J_CCCH_AL;
'CD_CE_HD2_HE2',     4, 2, 1, 3, J_HCCH_AL;       'CD_CE_HD2_HE3',     4, 2, 1, 3, J_HCCH_AL;
'CD_CE_HD3_HE2',     4, 2, 1, 3, J_HCCH_AL;       'CD_CE_HD3_HE3',     4, 2, 1, 3, J_HCCH_AL;
'CD_CG_HD2_HG2',     4, 2, 1, 3, J_HCCH_AL;       'CD_CG_HD2_HG3',     4, 2, 1, 3, J_HCCH_AL;
       
% Generic aromatics
'CE1_CE2_CZ_HE1',    4, 1, 3, 2, J_CCCH_AR;       'CE1_CE2_CZ_HE2',    4, 2, 3, 1, J_CCCH_AR;
'CE1_CZ_HE1_HZ',     4, 2, 1, 3, J_HCCH_AR;       'CE2_CZ_HE2_HZ',     4, 2, 1, 3, J_HCCH_AR
'CD2_CE2_CG_CZ',     3, 1, 2, 4, J_CCCC_AR;       'CB_CD1_CE1_CG',     1, 4, 2, 3, J_CCCC_AR;
'CD2_CE1_CE2_CZ',    1, 3, 4, 2, J_CCCC_AR;       'CD1_CE1_CG_CZ',     3, 1, 2, 4, J_CCCC_AR;
'CD1_CE1_CG_HE1',    4, 2, 1, 3, J_CCCH_AR;       'CD1_CE1_CZ_HD1',    4, 1, 2, 3, J_CCCH_AR;
'CD1_CE1_CZ_HZ',     4, 3, 2, 1, J_CCCH_AR;       'CD2_CE2_CG_HE2',    4, 2, 1, 3, J_CCCH_AR;
'CD2_CE2_CZ_HD2',    4, 1, 2, 3, J_CCCH_AR;       'CD2_CE2_CZ_HZ',     4, 3, 2, 1, J_CCCH_AR;
'CD2_CE2_HD2_HE2',   3, 1, 2, 4, J_HCCH_AR;       'CD1_CE1_HD1_HE1',   3, 1, 2, 4, J_HCCH_AR;
'CD1_CE1_CE2_CZ',    1, 2, 4, 3, J_CCCC_AR;       

% Generic amides
'CB_CG_HB2_ND2',     3, 1, 2, 4, [3.1 -0.6 0.4]; % http://dx.doi.org/10.1016/j.carres.2007.02.023
'CB_CG_HB3_ND2',     3, 1, 2, 4, [3.1 -0.6 0.4]; % http://dx.doi.org/10.1016/j.carres.2007.02.023
'CB_CG_HD21_ND2',    1, 2, 4, 3, [0.0 0.0 0.0];  % Tentative, needs checking
'CB_CG_HD22_ND2',    1, 2, 4, 3, [0.0 0.0 0.0];  % Tentative, needs checking
'CD_CG_HE21_NE2',    3, 4, 1, 2, [0.0 0.0 0.0];  % Tentative, needs checking
'CD_CG_HE22_NE2',    3, 4, 1, 2, [0.0 0.0 0.0];  % Tentative, needs checking
'CD_CG_HG2_NE2',     3, 2, 1, 4, [3.1 -0.6 0.4]; % http://dx.doi.org/10.1016/j.carres.2007.02.023
'CB_CD_CG_NE2',      1, 3, 2, 4, [0.0 0.0 0.0];  % Tentative, needs checking
'CA_CB_CG_ND2',      1, 2, 3, 4, [0.0 0.0 0.0];  % Tentative, needs checking
'CD_CG_HG3_NE2',     3, 2, 1, 4, [3.1 -0.6 0.4]; % http://dx.doi.org/10.1016/j.carres.2007.02.023 

% Extra parameters highlighted by the Boston visit
'CZ_HH22_NH1_NH2',   2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HH21_NH1_NH2',   2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HH12_NH1_NH2',   2, 3, 1, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HH11_NH1_NH2',   2, 3, 1, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HH22_NE_NH2',    2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HH21_NE_NH2',    2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HH12_NE_NH1',    2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HH11_NE_NH1',    2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HE_NE_NH2',      2, 3, 1, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CZ_HE_NE_NH1',      2, 3, 1, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CZ_NE_NH2',      1, 3, 2, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CZ_NE_NH1',      1, 3, 2, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CZ_HD3_NE',      3, 1, 4, 2, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CZ_HD2_NE',      3, 1, 4, 2, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_HD3_HE_NE',      2, 1, 4, 3, [4.22 -0.50  4.50]; % Tentative, needs checking
'CD_HD2_HE_NE',      2, 1, 4, 3, [4.22 -0.50  4.50]; % Tentative, needs checking
'CD_CG_CZ_NE',       2, 1, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CG_HE_NE',       2, 1, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CG_HG3_NE',      3, 2, 1, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CG_HG2_NE',      3, 2, 1, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CB_CD_CG_NE',       1, 3, 2, 4, [0.0 0.0 0.0]; % Tentative, needs checking

'CE_HE3_HZ3_NZ',     2, 1, 4, 3, [6.51 -1.76 1.60]; % Tentative, needs checking
'CE_HE2_HZ3_NZ',     2, 1, 4, 3, [6.51 -1.76 1.60]; % Tentative, needs checking
'CE_HE3_HZ2_NZ',     2, 1, 4, 3, [6.51 -1.76 1.60]; % Tentative, needs checking
'CE_HE2_HZ2_NZ',     2, 1, 4, 3, [6.51 -1.76 1.60]; % Tentative, needs checking
'CE_HE3_HZ1_NZ',     2, 1, 4, 3, [6.51 -1.76 1.60]; % Tentative, needs checking
'CE_HE2_HZ1_NZ',     2, 1, 4, 3, [6.51 -1.76 1.60]; % Tentative, needs checking
'CD_CE_HZ3_NZ',      1, 2, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CE_HZ2_NZ',      1, 2, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CE_HZ1_NZ',      1, 2, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CE_HD3_NZ',      3, 1, 2, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CE_HD2_NZ',      3, 1, 2, 4, [0.0 0.0 0.0]; % Tentative, needs checking
'CD_CE_CG_NZ',       3, 1, 2, 4, [0.0 0.0 0.0]; % Tentative, needs checking

'CA_H3_HA_N',        2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CA_H2_HA_N',        2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CA_H1_HA_N',        2, 4, 1, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CA_CB_H3_N',        2, 1, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CA_CB_H2_N',        2, 1, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'CA_CB_H1_N',        2, 1, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'C_CA_H3_N',         1, 2, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'C_CA_H2_N',         1, 2, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking
'C_CA_H1_N',         1, 2, 4, 3, [0.0 0.0 0.0]; % Tentative, needs checking

% Tryptophan side chain (GIAO DFT M06/cc-pVTZ data in PCM water)
'CH2_CZ3_HH2_HZ3'    3, 1, 2, 4, [0.0 0.0  9.41];   'CH2_CZ2_HH2_HZ2'    3, 1, 2, 4, [0.0 0.0 10.87];
'CH2_CZ2_CZ3_HZ3'    2, 1, 3, 4, [0.0 0.0  8.31];   'CH2_CZ2_CZ3_HZ2'    3, 1, 2, 4, [0.0 0.0  7.63];
'CE3_CZ3_HE3_HZ3'    3, 1, 2, 4, [0.0 0.0 11.09];   'CE3_CH2_CZ3_HH2'    1, 3, 2, 4, [0.0 0.0  8.12];
'CE3_CH2_CZ3_HE3'    2, 3, 1, 4, [0.0 0.0  8.56];   'CE3_CH2_CZ2_CZ3'    1, 4, 2, 3, [0.0 0.0  8.07];
'CE2_CH2_CZ2_HH2'    1, 3, 2, 4, [0.0 0.0  9.67];   'CE2_CH2_CZ2_CZ3'    1, 3, 2, 4, [0.0 0.0 10.44];
'CE2_CZ2_HZ2_NE1'    4, 1, 2, 3, [0.0 0.0  1.36];   'CE2_CZ2_HE1_NE1'    2, 1, 4, 3, [0.0 0.0  1.31];
'CE2_CH2_CZ2_NE1'    2, 3, 1, 4, [0.0 0.0  1.48];   'CD2_CE3_CZ3_HZ3'    1, 2, 3, 4, [0.0 0.0  8.32];
'CD2_CE3_CH2_CZ3'    1, 2, 4, 3, [0.0 0.0  8.85];   'CD2_CE2_CZ2_HZ2'    1, 2, 3, 4, [0.0 0.0  5.17];
'CD2_CE2_CH2_CZ2'    1, 2, 4, 3, [0.0 0.0  8.85];   'CD2_CE2_CE3_HE3'    2, 1, 3, 4, [0.0 0.0  7.92];
'CD2_CE2_CE3_CZ3'    2, 1, 3, 4, [0.0 0.0 10.44];   'CD2_CE2_CE3_CZ2'    3, 1, 2, 4, [0.0 0.0  8.07];
'CD2_CE2_HE1_NE1'    1, 2, 4, 3, [0.0 0.0  6.80];   'CD2_CE2_CE3_NE1'    3, 1, 2, 4, [0.0 0.0  1.01];
'CD1_HD1_HE1_NE1'    2, 1, 4, 3, [0.0 0.0  3.64];   'CD1_CE2_HD1_NE1'    2, 4, 1, 3, [0.0 0.0  6.87];
'CD1_CE2_CZ2_NE1'    3, 2, 4, 1, [0.0 0.0  2.34];   'CD1_CD2_CE2_NE1'    2, 3, 4, 1, [0.0 0.0  2.29];
'CD2_CE3_CG_HE3'     3, 1, 2, 4, [0.0 0.0  4.04];   'CD2_CE3_CG_CZ3'     3, 1, 2, 4, [0.0 0.0  4.69];
'CD2_CE2_CG_CZ2'     3, 1, 2, 4, [0.0 0.0  2.56];   'CD2_CE2_CG_NE1'     3, 1, 2, 4, [0.0 0.0  2.30];
'CD1_CG_HE1_NE1'     2, 1, 4, 3, [0.0 0.0  5.75];   'CD1_CE2_CG_NE1'     3, 1, 4, 2, [0.0 0.0  3.34];
'CD1_CD2_CE3_CG'     1, 4, 2, 3, [0.0 0.0  5.66];   'CD1_CD2_CG_NE1'     2, 3, 1, 4, [0.0 0.0  3.54];
'CB_CD2_CE3_CG'      1, 4, 2, 3, [0.0 0.0  1.70];   'CB_CD1_CG_NE1'      1, 3, 2, 4, [0.0 0.0  1.21];

};      

% Collision detection exceptions (TRP, HIS and PRO, where two-bond and three-bond couplings collide)
allowed_collision_list={'CD2_CE1_ND1_NE2','CE1_CG_ND1_NE2','CD2_CE1_CG_ND1','CD1_CE1_CG_CZ', 'CD2_CE2_CE3_CZ2',...
                        'CD1_CD2_CE2_CG','CA_CB_CD_CG','CA_CD_CG_N','CA_CB_CD_N','CD2_CE1_CG_NE2', 'CD2_CE2_CG_NE1',...
                        'CD2_CE2_CH2_CZ2','CD2_CE2_CE3_CZ3','CD1_CE2_CG_NE1','CA_CB_CG2_HA', 'CD1_CD2_CG_NE1',...
                        'CA_CB_CG1_HA','CD2_CE2_CG_CZ','CD1_CE1_CE2_CZ','CD2_CE3_CH2_CZ3','CE2_CH2_CZ2_CZ3',...
                        'CE3_CH2_CZ2_CZ3'};

% Loop over subgraphs
disp(' '); disp('############# THREE-BOND J-COUPLING SUMMARY ###############');
for n=1:size(subgraphs,1)
    
    % Extract labels, residues and numbers
    spin_labels=pdb_id(subgraphs(n,:));
    spin_numbers=numbers(subgraphs(n,:));
    spin_resnames=aa_typ(subgraphs(n,:));
    spin_resnums=aa_num(subgraphs(n,:));
    
    % Index the spins
    [spin_labels,index]=sort(spin_labels);
    spin_numbers=spin_numbers(index);
    spin_resnames=spin_resnames(index);
    spin_resnums=spin_resnums(index);
    
    % Consult the database
    if isscalar(spin_labels)
        
        % Inform the user if no J-coupling found
        disp(['WARNING: solitary spin detected, ' spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1}]);
        
    elseif numel(spin_labels)==2
        
        % Inform the user if no J-coupling found
        disp(['WARNING: isolated spin pair detected, ' spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1} ' and ' ...
                                                       spin_resnames{2} '(' num2str(spin_resnums(2)) '):' spin_labels{2}]);
    
    elseif numel(spin_labels)==3
        
        % Inform the user if no J-coupling found
        disp(['WARNING: isolated spin triad detected, ' spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1} ', '...
                                                        spin_resnames{2} '(' num2str(spin_resnums(2)) '):' spin_labels{2} ' and '...
                                                        spin_resnames{3} '(' num2str(spin_resnums(3)) '):' spin_labels{3}]);
    else
    
        % Assemble specifiction string
        spec_strings=quads_database(strcmp([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4}],quads_database(:,1)),:);
        
        % Complain if nothing found
        if isempty(spec_strings)
            
            disp(['WARNING: unknown atom quad or unphysical proximity, ' spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4} ' in amino acids '...
                                                                         spin_resnames{1} '(' num2str(spin_resnums(1)) '), '...
                                                                         spin_resnames{2} '(' num2str(spin_resnums(2)) '), '...
                                                                         spin_resnames{3} '(' num2str(spin_resnums(3)) '), '...
                                                                         spin_resnames{4} '(' num2str(spin_resnums(4)) ')']);
        else
            
            % Loop over specifications
            for k=1:size(spec_strings,1)
                
                % Extract spin numbers
                spin_a=spin_numbers(spec_strings{k,2}); spin_b=spin_numbers(spec_strings{k,3});
                spin_c=spin_numbers(spec_strings{k,4}); spin_d=spin_numbers(spec_strings{k,5});
                
                % Match bonding pattern
                if proxmatrix(spin_a,spin_b)&&proxmatrix(spin_b,spin_c)&&proxmatrix(spin_c,spin_d)
                    
                    % Get Karplus coefficients
                    A=spec_strings{k,6}(1); B=spec_strings{k,6}(2); C=spec_strings{k,6}(3);
                    
                    % Get dihedral angle
                    theta=dihedral(coords{spin_a},coords{spin_b},coords{spin_c},coords{spin_d});
                    
                    % Assign the J-coupling
                    if (~isempty(jmatrix{spin_a,spin_d}))&&(~ismember([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4}],allowed_collision_list))
                        
                        % Complain and bomb out
                        error(['atom quad: ' spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4} ' in amino acids '...
                                             spin_resnames{1} '(' num2str(spin_resnums(1)) '), '...
                                             spin_resnames{2} '(' num2str(spin_resnums(2)) '), '...
                                             spin_resnames{3} '(' num2str(spin_resnums(3)) '), '...
                                             spin_resnames{4} '(' num2str(spin_resnums(4)) ') collides with other data.']);
                        
                    else
                        
                        % Write the array
                        jmatrix{spin_a,spin_d}=A*cosd(theta)^2+B*cosd(theta)+C;
                        
                        % Inform the user
                        disp(['Estimated three-bond J-coupling from ' pad([aa_typ{spin_a} '(' num2str(aa_num(spin_a)) '):' pdb_id{spin_a}],19)   ' to '...
                                                                      pad([aa_typ{spin_d} '(' num2str(aa_num(spin_d)) '):' pdb_id{spin_d}],19) ' '...
                                                                      pad(num2str(jmatrix{spin_a,spin_d},'%5.1f'),6,'left') ' Hz']);
                                                                  
                    end
                    
                end
                
            end
            
        end
        
    end

end

end

% Consistency enforcement
function grumble(aa_num,aa_typ,pdb_id,coords)
if ~isnumeric(aa_num)
    error('aa_num must be a vector of positive integers.');
end
if ~iscell(aa_typ)
    error('aa_typ must be a cell array of strings.');
end
if ~iscell(pdb_id)
    error('pdb_id must be a cell array of strings.');
end
if ~iscell(coords)
    error('coords must be a cell array of 3-vectors.');
end
if (numel(aa_num)~=numel(aa_typ))||...
   (numel(aa_typ)~=numel(pdb_id))||...
   (numel(pdb_id)~=numel(coords))
    error('all four inputs must have the same number of elements.');
end
end

% Niels Bohr, on Dirac's equation:
%
% "Simply put an explanation of the theory on a poster, tack it up on a tree
% in the jungle, and any elephant (a beast noted for its wisdom) that passed
% by would immediately become so engrossed trying to figure it out that it
% could be packed up and delivered to the Copenhagen Zoo before it realized
% anything had happened."

