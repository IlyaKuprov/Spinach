% RNA assignments of J-couplings from literature values and Karplus cur-
% ves. Syntax:
%
%           jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)
%
% Parameters:
%
%       nuc_num   - a vector of nucleotide numbers
%
%       nuc_typ   - a cell array of nucleotide types
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
% i.kuprov@soton.ac.uk
% ztw1e12@soton.ac.uk (Zenawi Welderufael)
%
% <https://spindynamics.org/wiki/index.php?title=Guess_j_nuc.m>

function jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)

% Check consistency
grumble(nuc_num,nuc_typ,pdb_id,coords);

% Preallocate the answer
jmatrix=cell(numel(coords),numel(coords));

% Number the atoms
numbers=1:numel(coords);

% Determine the connectivity graph
proxmatrix=false(numel(coords),numel(coords));
for n=1:numel(coords)
    for k=1:numel(coords)
        if (norm(coords{n}-coords{k},2)<1.55), proxmatrix(n,k)=1; end
    end
end

% Get all connected subgraphs up to size 2
subgraphs=dfpt(sparse(proxmatrix),2);

% Set generic coupling values
J_NH=86.0; J_CN=20.0; J_CH=180.0; J_CC=65.0;

% Spec all ordered connected pairs
pairs_database={...

% RNA specific values from the literature
'C1p_N9'   , 11.0;    'C1p_N1'   , 12.0;    'H_N1'     , 95.0;   'H3_N3'    , 91.0;   
'C6_N1'    , 7.5;     'C4_N3'    , 10.7;    'C4_C5'    , 65.0;   'C5_H5'    , 176;  
'C8_H8'    , 216.0;   'C6_H6'    , 185.0;   'C5_C6'    , 67.0;   'C2_N3'    , 19.0;  
'C2_N1'    , 19.0;    'C4_N9'    , 20.0;    'H41_N4'   , 86.0;   'C4_N4'    , 20.0;
'C8_N9'    , 11.0;    'C1p_H1p'  , 170.0;   'C_C5'     , 55.0;   'H42_N4'   , 86.0;    

% Generics for the conjugated ring (to be replaced with DFT values - Zenawi)
'C2_H2'    , J_CH;    'C5_N7'    , J_CN;    'C8_N7'    , J_CN;   'C6_N6'    , J_CN;
'H61_N6'   , J_NH;    'H62_N6'   , J_NH;    'C2_N2'    , J_CN;   'H21_N2'   , J_NH;
'H1_N1'    , J_NH;    'H22_N2'   , J_NH;    'H61_N1'   , J_NH;

% Generics for the sugar ring (to be replaced with DFT values - Zenawi)
'C1p_C2p'  , J_CC;    'C2p_H2p'  , J_CH;    'C3p_C4p'  , J_CC;   'C4p_H4p'  , J_CH;
'C4p_C5p'  , J_CC;    'C5p_H5p'  , J_CH;    'C5p_H5pp' , J_CH;   'C2p_C3p'  , J_CC;   
'C3p_H3p'  , J_CH;    'C2p_H2pp' , J_CH;     'C2p_H2p1' , J_CH;  'C5p_H5p2' , J_CH
'C5p_H5p1' , J_CH};      

% Loop over subgraphs
disp(' '); disp('############# ONE-BOND J-COUPLING SUMMARY ###############');

for n=1:size(subgraphs,1)
    
    % Extract labels, residues and numbers
    spin_labels=pdb_id(subgraphs(n,:));
    spin_numbers=numbers(subgraphs(n,:));
    spin_resnames=nuc_typ(subgraphs(n,:));
    spin_resnums=nuc_num(subgraphs(n,:));
    
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
J_NCN=0;     % Literature reference
J_NCH=0;     % Literature reference
J_CNH=0;     % Literature reference
J_CNC=0;     % Literature reference
J_HCH=-12.0; % Typical for CH2
J_HNH=-10.0; % Typical for NH2

% Spec all ordered connected triples (sorted atoms, bonding order, coupling from atom 1 to atom 3)
triples_database={...
    
% Backbone data from the literature
'C4p_C5p_H5p'  1, 2, 3,  -5.5;   'C4p_C5p_H5pp'  1, 2, 3,   2.0;   'C4p_C5p_H4p' 3, 1, 2,  -5.5;
'C2p_C3p_H3p'  1, 2, 3,  -2.3;   'C2p_C3p_H2p'   3, 2, 1,   2.3;  

% Conjugated ring data from the literature
'C4_C5_C6'     1, 2, 3,   9.5;   'C4_C8_N9'      1, 3, 2,   8.0;   'C2_C4_N3'    1, 3, 2,  10.0;
'C8_H8_N9'     2, 1, 3,   8.0;   'C8_H8_N7'      2, 1, 3,  11.0;   'C2_H2_N1'    2, 1, 3,  15.0;
'C2_H2_N3'     2, 1, 3,  15.0;

% Generic numbers for the conjugated ring (to be replaced with DFT values - Zenawi)
'C6_N1_N6'     2, 1, 3,  J_NCN;  'C2_N1_N3'      2, 1, 3,  J_NCN;  'C2_C6_N1'     1, 3, 2,  J_CNC;                                     
'C4_C5_N3'     2, 1, 3,  J_CCN;  'C4_N3_N9'      2, 1, 3,  J_NCN;  'C1p_C4_N9'    1, 3, 2,  J_CNC;  
'C4_C5_N9'     2, 1, 3,  J_CCN;  'C5_C6_N6'      1, 2, 3,  J_CCN;  'C5_C8_N7'     1, 3, 2,  J_CNC;   
'C4_C5_N7'     1, 2, 3,  J_CCN;  'C5_C6_N7'      2, 1, 3,  J_CCN;  'C6_H61_N6'    1, 3, 2,  J_CNH;  
'C6_H62_N6'    1, 3, 2,  J_CNH;  'C8_N7_N9'      2, 1, 3,  J_NCN;  'C1p_C8_N9'    2, 3, 1,  J_CNC;  
'C2_N1_N2'     2, 1, 3,  J_NCN;  'C2_H1_N1'      1, 3, 2,  J_CNH;  'C2_H21_N2'    1, 3, 2,  J_CNH;  
'C2_H22_N2'    1, 3, 2,  J_CNH;  'C2_N2_N3'      3, 1, 2,  J_NCN;  'C6_H1_N1'     1, 3, 2,  J_CNH;      
'C5_C6_N1'     2, 1, 3,  J_CCN;  'C6_H6_N1'      3, 1, 2,  J_NCH;  'C1p_C2_N1'    2, 3, 1,  J_CNC;                                    
'C4_N3_N4'     2, 1, 3,  J_NCN;  'C4_H41_N4'     1, 3, 2,  J_CNH;  'C4_H42_N4'    1, 3, 2,  J_CNH;  
'C4_C5_N4'     2, 1, 3,  J_CCN;  'C5_C6_H6'      1, 2, 3,  J_CCH;  'C1p_C6_N1'    2, 3, 1,  J_CNC;
'C5_C6_H5'     2, 1, 3,  J_CCH;  'C2_H3_N3'      1, 3, 2,  J_CNH;  'C4_H3_N3'     1, 3, 2,  J_CNH;    
'C4_C5_H5'     1, 2, 3,  J_CCH;  'C1p_H1p_N9'    2, 1, 3,  J_NCH;  'C1p_C2p_N9'   2, 1, 3,  J_CCN;
'H21_H22_N2'   1, 3, 2,  J_HNH;  'C1p_H1p_N1'    2, 1, 3,  J_NCH;  'C1p_C2p_N1'   2, 1, 3,  J_CCN;
'H41_H42_N4'   1, 3, 2,  J_HNH;  'H61_H62_N6'    1, 3, 2,  J_HNH; 

% Generic numbers for the backbone (to be replaced with DFT values - Zenawi)
'C1p_C2p_C3p'   1, 2, 3,  J_CCC;  'C1p_C2p_H2p'   1, 2, 3,  J_CCH;   'C2p_C3p_C4p'  1, 2, 3,  J_CCC;
'C3p_C4p_H4p'   3, 1, 2,  J_CCH;  'C3p_C4p_C5p'   1, 2, 3,  J_CCC;   'C2p_C3p_H2p'  1, 2, 3,  J_CCH;  
'C3p_C4p_H3p'   2, 1, 3,  J_CCH;  'C5p_H5p_H5pp'  2, 1, 3,  J_HCH;   'C1p_C2p_H1p'  3, 1, 2,  J_CCH;
'C1p_C2p_H2pp'  3, 2, 1,  J_CCH;  'C2p_C3p_H2pp'  3, 1, 2,  J_CCH;  

'C2p_C3p_H2p1' 1, 2, 3,  J_CCH; 'C1p_C2p_H2p1'  1, 2, 3,  J_CCH;  'C4p_C5p_H5p2'  1, 2, 3,  J_CCH;  
'C4p_C5p_H5p1' 1, 2, 3, J_CCH; 'C5p_H5p1_H5p2'  2, 1, 3, J_HCH};

% Exception list for four-membered rings
allowed_collision_list={};

% Loop over subgraphs
disp(' '); disp('############# TWO-BOND J-COUPLING SUMMARY ###############');
for n=1:size(subgraphs,1)
    
    % Extract labels, residues and numbers
    spin_labels=pdb_id(subgraphs(n,:));
    spin_numbers=numbers(subgraphs(n,:));
    spin_resnames=nuc_typ(subgraphs(n,:));
    spin_resnums=nuc_num(subgraphs(n,:));
    
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
            disp(['WARNING: unknown atom triple or unphysical proximity, ' spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} ' in residues '...
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
                disp(['Estimated two-bond J-coupling from ' pad([nuc_typ{spin_a} '(' num2str(nuc_num(spin_a)) '):' pdb_id{spin_a}],19) ' to '...
                                                            pad([nuc_typ{spin_c} '(' num2str(nuc_num(spin_c)) '):' pdb_id{spin_c}],19) ' '...
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

% List of ordered connected quads (sorted atoms, bonding order, Karplus coefficients)
quads_database={...
'C3p_C4p_C5p_H5p',  4, 3, 2, 1, [0.0 0.0 0.0];  'C4p_C5p_H4p_H5p',   4, 2, 1, 3, [0.0 0.0 10.42];
'C3p_C4p_C5p_H5pp', 4, 3, 2, 1, [0.0 0.0 0.0];  'C4p_C5p_H4p_H5pp',  4, 2, 1, 3, [0.0 0.0 8.51];
'C1p_C8_H8_N9',     1, 4, 2, 3, [0.0 0.0 0.0];  'C1p_C8_N7_N9',      1, 4, 2, 3, [0.0 0.0 0.0];
'C4_C8_H8_N9',      1, 4, 2, 3, [0.0 0.0 0.0];  'C5_C8_N7_N9',       4, 2, 3, 1, [0.0 0.0 0.0];  
'C4_C5_C8_N7',      3, 4, 2, 1, [0.0 0.0 0.0];  'C5_C6_C8_N7',       3, 4, 1, 2, [0.0 0.0 9.00];
'C5_C6_N6_N7',      4, 1, 2, 3, [0.0 0.0 0.0];  'C5_C6_N1_N7',       4, 1, 2, 3, [0.0 0.0 0.0];
'C4_C5_C6_N6',      1, 2, 3, 4, [0.0 0.0 0.0];  'C4_C5_C6_N1',       1, 2, 3, 4, [0.0 0.0 0.0];
'C5_C6_H61_N6',     1, 2, 4, 3, [0.0 0.0 0.0];  'C5_C6_H62_N6',      1, 2, 4, 3, [0.0 0.0 0.0];
'C6_H61_N1_N6',     3, 1, 4, 2, [0.0 0.0 0.0];  'C6_H62_N1_N6',      3, 1, 4, 2, [0.0 0.0 0.0];
'C2_C6_N1_N6',      1, 3, 2, 4, [0.0 0.0 0.0];  'C2_C6_N1_N3',       4, 1, 3, 2, [0.0 0.0 0.0];  
'C2_C4_H2_N3',      2, 4, 1, 3, [0.0 0.0 0.0];  'C2_C6_H2_N1',       3, 1, 4, 2, [0.0 0.0 0.0];
'C2_C4_C5_N3',      3, 2, 4, 1, [0.0 0.0 0.0];  'C2_C4_N3_N9',       4, 2, 3, 1, [0.0 0.0 0.0];
'C4_C5_C6_N3',      3, 2, 1, 4, [0.0 0.0 0.0];  'C4_C5_N3_N7',       4, 2, 1, 3, [0.0 0.0 0.0];
'C4_C5_C6_N9',      3, 2, 1, 4, [0.0 0.0 0.0];  'C4_C8_N3_N9',       3, 1, 4, 2, [0.0 0.0 0.0];  
'C2_C5_C6_N1',      1, 4, 3, 2, [0.0 0.0 0.0];  'C5_C6_H1_N1',       3, 4, 2, 1, [0.0 0.0 0.0];
'C2_C6_N1_N2',      4, 1, 3, 2, [0.0 0.0 0.0];  'C2_H1_N1_N2',       4, 1, 3, 2, [0.0 0.0 0.0];
'C2_H1_N1_N3',      4, 1, 3, 2, [0.0 0.0 0.0];  'C2_H21_N2_N3',      4, 1, 3, 2, [0.0 0.0 0.0];  
'C2_H21_N1_N2',     3, 1, 4, 2, [0.0 0.0 0.0];  'C2_H22_N2_N1',      3, 1, 4, 2, [0.0 0.0 0.0];
'C2_C4_N2_N3',      2, 4, 1, 3, [0.0 0.0 0.0];  'C1p_C4_C5_N9',      3, 2, 4, 1, [0.0 0.0 0.0];
'C1p_C4_N3_N9',     3, 2, 4, 1, [0.0 0.0 0.0];  'C1p_C5_C6_N1',      1, 4, 3, 2, [0.0 0.0 0.0];
'C2_C6_H6_N1',      1, 4, 2, 3, [0.0 0.0 0.0];  'C5_C6_H5_N1',       4, 2, 1, 3, [0.0 0.0 0.0];
'C4_C5_C6_H6',      4, 3, 2, 1, [0.0 0.0 0.0];  'C5_C6_H5_H6',       4, 2, 1, 3, [0.0 0.0 7.00];
'C4_C5_C6_N4',      3, 2, 1, 4, [0.0 0.0 0.0];  'C4_C5_H5_N3',       3, 2, 1, 4, [0.0 0.0 0.0];  
'C4_C5_H5_N4',      3, 2, 1, 4, [0.0 0.0 0.0];  'C2_H22_N2_N3',      4, 1, 3, 2, [0.0 0.0 0.0];
'C4_C5_H41_N4',     2, 1, 4, 3, [0.0 0.0 0.0];  'C4_C5_H42_N4',      2, 1, 4, 3, [0.0 0.0 0.0];
'C4_H41_N3_N4',     3, 1, 4, 2, [0.0 0.0 0.0];  'C4_H42_N3_N4',      3, 1, 4, 2, [0.0 0.0 0.0];
'C2_C4_N3_N4',      4, 2, 3, 1, [0.0 0.0 0.0];  'C1p_C2_N1_N3',      4, 2, 3, 1, [0.0 0.0 0.0];  
'C1p_C6_H6_N1',     1, 4, 2, 3, [0.0 0.0 0.0];  'C4_C5_H3_N3',       2, 1, 4, 3, [0.0 0.0 0.0];  
'C2_H3_N1_N3',      2, 4, 1, 3, [0.0 0.0 0.0];  'C5_C8_H8_N7',       3, 2, 4, 1, [0.0 0.0 0.0];
'C1p_C2p_H2p_N1',   3, 2, 1, 4, [0.0 0.0 0.0];  'C1p_C2p_H2p_N9',    3, 2, 1, 4, [0.0 0.0 0.0];
'C1p_C2p_C3p_N1',   3, 2, 1, 4, [0.0 0.0 0.0];  'C1p_C2p_C3p_N9',    3, 2, 1, 4, [0.0 0.0 0.0]; 
'C2p_C3p_C4p_C5p',  4, 3, 2, 1, [0.0 0.0 0.0];  'C3p_C4p_H3p_H4p',   4, 2, 1, 3, [0.0 0.0 5.79];
'C2p_C3p_C4p_H4p',  4, 3, 2, 1, [0.0 0.0 0.0];  'C3p_C4p_C5p_H3p',   3, 2, 1, 4, [0.0 0.0 0.0]; 
'C1p_C2p_C3p_C4p',  4, 3, 2, 1, [0.0 0.0 0.0];  'C2p_C3p_H2p_H3p',   4, 2, 1, 3, [0.0 0.0 0.0];
'C1p_C2p_C3p_H3p',  4, 3, 2, 1, [0.0 0.0 0.0];  'C2p_C3p_C4p_H2p',   3, 2, 1, 4, [0.0 0.0 0.0]; 
'C2p_C3p_H2pp_H3p', 3, 1, 2, 4, [0.0 0.0 5.82]; 'C2_C4_N1_N3',       3, 1, 4, 2, [0.0 0.0 0.0];
'C1p_C2p_H2pp_N1',  4, 1, 2, 3, [0.0 0.0 0.0];  'C1p_C2_C2p_N1',     2, 4, 1, 3, [0.0 0.0 0.0];
'C1p_C2p_C6_N1',    3, 4, 1, 2, [0.0 0.0 0.0];  'C1p_C2p_C3p_H1p',   4, 1, 2, 3, [0.0 0.0 0.0];
'C1p_C2p_H1p_H2pp', 3, 1, 2, 4, [0.0 0.0 4.05]; 'C1p_C2_H1p_N1',     2, 4, 1, 3, [0.0 0.0 0.0];
'C1p_C6_H1p_N1',    2, 4, 1, 3, [0.0 0.0 0.0];  'C2p_C3p_C4p_H2pp',  4, 1, 2, 3, [0.0 0.0 0.0];
'C2_H22_N1_N2',     3, 1, 4, 2, [0.0 0.0 0.0];  'C4_C5_C8_N9',       3, 4, 1, 2, [0.0 0.0 0.0];
'C4_C5_N7_N9',      4, 1, 2, 3, [0.0 0.0 0.0];  'C1p_C2p_H2pp_N9',   4, 1, 2, 3, [0.0 0.0 0.0];
'C1p_C2p_C8_N9',    3, 4, 1, 2, [0.0 0.0 0.0];  'C1p_C2p_C4_N9',     3, 4, 1, 2, [0.0 0.0 0.0];
'C1p_C8_H1p_N9',    2, 4, 1, 3, [0.0 0.0 0.0];  'C1p_C4_H1p_N9',     2, 4, 1, 3, [0.0 0.0 0.0];
'C4_C8_N7_N9',      3, 2, 4, 1, [0.0 0.0 0.0];  'C2p_C3p_H2p1_H3p',  4, 2, 1, 3, [0.0 0.0 0.0];
'C1p_C2p_H2p1_N1',  3, 2, 1, 4, [0.0 0.0 0.0];  'C1p_C2p_H1p_H2p1'   3, 1, 2, 4, [0.0 0.0 4.05];
'C2p_C3p_C4p_H2p1', 3, 2, 1, 4, [0.0 0.0 0.0];  'C3p_C4p_C5p_H5p2',  4, 3, 2, 1, [0.0 0.0 0.0];
'C3p_C4p_C5p_H5p1', 4, 3, 2, 1, [0.0 0.0 0.0];  'C4p_C5p_H4p_H5p2',  4, 2, 1, 3, [0.0 0.0 10.42];
'C1p_C2p_H2p1_N9',  3, 2, 1, 4, [0.0 0.0 0.0];  'C4p_C5p_H4p_H5p1',  4, 2, 1, 3, [0.0 0.0 10.42];
};

% Exception list for five-membered rings
allowed_collision_list={'C4_C5_C6_N1','C4_C5_N7_N9','C1p_C4_H1p_N9'};

% Loop over subgraphs
disp(' '); disp('############# THREE-BOND J-COUPLING SUMMARY ###############');
for n=1:size(subgraphs,1)
    
    % Extract labels, residues and numbers
    spin_labels=pdb_id(subgraphs(n,:));
    spin_numbers=numbers(subgraphs(n,:));
    spin_resnames=nuc_typ(subgraphs(n,:));
    spin_resnums=nuc_num(subgraphs(n,:));
    
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
        spec_string=quads_database(strcmp([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4}],quads_database(:,1)),:);
        
        % Complain if nothing found
        if isempty(spec_string)
            
            disp(['WARNING: unknown atom quad or unphysical proximity, ' spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4} ' in residues '...
                                                                         spin_resnames{1} '(' num2str(spin_resnums(1)) '), '...
                                                                         spin_resnames{2} '(' num2str(spin_resnums(2)) '), '...
                                                                         spin_resnames{3} '(' num2str(spin_resnums(3)) '), '...
                                                                         spin_resnames{4} '(' num2str(spin_resnums(4)) ')']);
        elseif size(spec_string,1)>1
            
            % Complain if multiple records found
            error([spec_string{1} ' coupling has been specified multiple times in the coupling database.']);
            
        else
                
            % Extract spin numbers
            spin_a=spin_numbers(spec_string{2}); spin_b=spin_numbers(spec_string{3});
            spin_c=spin_numbers(spec_string{4}); spin_d=spin_numbers(spec_string{5});
            
            % Get Karplus coefficients
            A=spec_string{6}(1); B=spec_string{6}(2); C=spec_string{6}(3);
            
            % Get dihedral angle
            theta=dihedral(coords{spin_a},coords{spin_b},coords{spin_c},coords{spin_d});
                    
            % Assign the J-coupling
            if (~isempty(jmatrix{spin_a,spin_d}))&&(~ismember([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4}],allowed_collision_list))
                
                % Complain and bomb out
                error(['atom quad: ' spin_labels{1} '_' spin_labels{2} '_' spin_labels{3} '_' spin_labels{4} ' in residues '...
                       spin_resnames{1} '(' num2str(spin_resnums(1)) '), '...
                       spin_resnames{2} '(' num2str(spin_resnums(2)) '), '...
                       spin_resnames{3} '(' num2str(spin_resnums(3)) '), '...
                       spin_resnames{4} '(' num2str(spin_resnums(4)) ') collides with other data.']);
                
            else
                
                % Write the array
                jmatrix{spin_a,spin_d}=A*cosd(theta)^2+B*cosd(theta)+C;
                
                % Inform the user
                disp(['Estimated three-bond J-coupling from ' pad([nuc_typ{spin_a} '(' num2str(nuc_num(spin_a)) '):' pdb_id{spin_a}],19)   ' to '...
                      pad([nuc_typ{spin_d} '(' num2str(nuc_num(spin_d)) '):' pdb_id{spin_d}],19) ' ' pad(num2str(jmatrix{spin_a,spin_d},'%5.1f'),6,'left') ' Hz']);
                
            end
            
        end
                
    end

end

end

% Consistency enforcement
function grumble(nuc_num,nuc_typ,pdb_id,coords)
if ~isnumeric(nuc_num)
    error('nuc_num must be a vector of positive integers.');
end
if ~iscell(nuc_typ)
    error('nuc_typ must be a cell array of strings.');
end
if ~iscell(pdb_id)
    error('pdb_id must be a cell array of strings.');
end
if ~iscell(coords)
    error('coords must be a cell array of 3-vectors.');
end
if (numel(nuc_num)~=numel(nuc_typ))||...
   (numel(nuc_typ)~=numel(pdb_id))||...
   (numel(pdb_id)~=numel(coords))
    error('all four inputs must have the same number of elements.');
end
end

% It's a terrible thing, I think, in life to wait until you are 
% ready. I have this feeling that actually no one is ever ready
% to do anything. There is almost no such thing as ready. There
% is only now.
%
% Hugh Laurie

