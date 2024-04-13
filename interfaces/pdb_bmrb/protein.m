% Protein data import function. Parses PDB and BMRB data, runs a J-coupl-
% ing guess, a CSA guess and outputs Spinach data structures. Syntax:
%
%            [sys,inter]=protein(pdb_file,bmrb_file,options)
%
% Parameters:
%
%        pdb_file   - string containing the name of the PDB file
%
%       bmrb_file   - string containing the name of the BMRB file
%
% options.select    - 'backbone' imports protein backbone up to
%                     CB and HB, 'backbone-minimal' only imports
%                     the backbone, 'backbone-hsqc' is the same
%                     as backbone, but with GLN and ASN side chain
%                     amide groups included, 'all' imports every-
%                     thing that is assigned in BMRB. If a list of
%                     numbers is supplied, spins with those num-
%                     bers in the PDB file are imported, but only
%                     if they are assigned in the PDB.
%
% options.pdb_mol   - the number of molecule if there are multiple 
%                     molecules in the pdb file 
%
% options.noshift   - 'keep' places unassigned atoms between -1 and
%                     0 ppm, 'delete' removes them from the system
%
% options.deuterate - a cell array of character strings, replaces 
%                     protons with the specified PDB identifiers 
%                     with deuterons
%
% Outputs:
%
%    sys.isotopes          - Nspins x 1 cell array of strings
% 
%    sys.labels            - Nspins x 1 cell array of strings containing 
%                            standard IUPAC protein atom labels
% 
%    inter.coordinates     - Nspins x 3 matrix, Angstrom.
% 
%    inter.zeeman.iso      - Nspins x 1 cell array of numbers, ppm. 
%                            Isotropic chemical shifts go here.
% 
%    inter.zeeman.matrix   - Nspins x 1 cell array of 3x3 matrices, ppm. 
%                            Chemical shift anisotropies go here.
% 
%    inter.coupling.scalar - Nspins x Nspins cell array of scalar coup-
%                            lings, all in Hz.
%
%    aux.pdb_aa_num        - pdb amino acid number for each spin
%
%    aux.pdb_aa_typ        - pdb amino acid type for each spin
%
% i.kuprov@soton.ac.uk
% zenawi.welderufael@soton.ac.uk
% andras_boeszoermenyi@hms.harvard.edu
% matt.walker@soton.ac.uk
% m.g.concilio@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=protein.m>

function [sys,inter,aux]=protein(pdb_file,bmrb_file,options)

% Set defaults
if ~isfield(options,'deuterate'), options.deuterate={}; end

% Check consistency
grumble(pdb_file,bmrb_file,options);

% Parse the PDB file
[pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords]=read_pdb_pro(pdb_file,options.pdb_mol);

% Parse the BMRB file
[bmrb_aa_num,bmrb_aa_typ,bmrb_atom_id,bmrb_chemsh]=read_bmrb(bmrb_file);

% Remove oxygens, sulphurs and terminal atoms
kill_mask=ismember(pdb_atom_id,{'O','OE','OE1','OE2','OD1','OD2','OG','OG1','HG1',...
                                'OG2','OH','HH','SD','SG','OXT','O''','O'''''});
pdb_aa_num(kill_mask)=[]; pdb_atom_id(kill_mask)=[]; 
pdb_aa_typ(kill_mask)=[]; pdb_coords(kill_mask)=[];
disp('WARNING: oxygen, sulphur, and OH protons will not appear in the simulation.');

% Match chemical shifts
pdb_chemsh=cell(numel(pdb_atom_id),1);
disp(' '); disp('############# CHEMICAL SHIFT IMPORT SUMMARY ###############');
for n=1:numel(pdb_atom_id)
    
    % Pull the current amino acid from BMRB
    select_mask=(bmrb_aa_num==pdb_aa_num(n));
    bmrb_atoms=bmrb_atom_id(select_mask);
    bmrb_shifts=bmrb_chemsh(select_mask);
    
    % Make sure amino acid types match
    if ~all(strcmp(pdb_aa_typ{n},bmrb_aa_typ(select_mask)))
        disp(['Amino acid residue ' num2str(n)]);
        disp(['PDB file says  ' pdb_aa_typ{n}]);
        disp('BMRB file says:'); disp(bmrb_aa_typ(select_mask));
        error('Amino acid type mismatch between PDB and BMRB files.');
    end
    
    % Ugly heuristics (sorry!) to match PDB and BMRB data
    if ismember(pdb_atom_id{n},bmrb_atoms)
        
        % If the atom is found, just import the shift
        pdb_chemsh{n}=bmrb_shifts(strcmp(pdb_atom_id{n},bmrb_atoms));
        disp(['Imported chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right')...
                                             pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                         
    elseif strcmp(pdb_atom_id{n},'N')&&(pdb_aa_num(n)==1)
        
        % Ignore the N-terminal amino group nitrogen (usually hard to assign)
        disp(['WARNING: N-terminal NH2 group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
        
    elseif ismember(pdb_atom_id{n},{'HG'})&&ismember(pdb_aa_typ(n),{'SER'})
        
        % Ignore OH protons in serine (usually rapid exchange)
        disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
        
    elseif ismember(pdb_atom_id{n},{'HG1'})&&ismember(pdb_aa_typ(n),{'THR'})
        
        % Ignore OH protons in threonine (usually rapid exchange)
        disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
        
    elseif ismember(pdb_atom_id{n},{'HB3'})&&ismember('HB2',bmrb_atoms)&&...
           ismember(pdb_aa_typ{n},{'MET','GLU','GLN','LYS','LEU','SER','HIS','ARG'})
        
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HB2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
        
    elseif ismember(pdb_atom_id{n},{'HB1','HB2','HB3'})&&...
           ismember('HB',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ALA'})
       
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HB',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HG21','HG22','HG23','HG1','HG3'})&&...
           ismember('HG2',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ILE','THR','VAL','GLU','LYS','PRO','GLN','ARG'})
        
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HG2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HG11','HG12','HG13','HG2','HG3'})&&...
           ismember('HG1',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'VAL'})
       
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HG1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HG13'})&&...
           ismember('HG12',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ILE'})
       
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HG12',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HD11','HD12','HD13','HD2','HD3'})&&...
           ismember('HD1',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ILE','LEU'})
       
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HD1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HD21','HD22','HD23','HD1','HD3'})&&...
           ismember('HD2',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'LYS','LEU','PRO','ARG'})
       
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HD2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HE3'})&&...
           ismember('HE2',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'LYS'})
        
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HE2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HE1','HE2','HE3'})&&...
           ismember('HE',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'MET'})
       
        % When chemical shift is given for just one proton of a methyl
        % group, use it for all three protons (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HE',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'CD2'})&&ismember('CD1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        
        % When chemical shift is given for one proton in a symmetric
        % aromatic ring, use it for the other one (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('CD1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'CE2'})&&ismember('CE1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        
        % When chemical shift is given for one proton in a symmetric
        % aromatic ring, use it for the other one (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('CE1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HD2'})&&ismember('HD1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        
        % When chemical shift is given for one proton in a symmetric
        % aromatic ring, use it for the other one (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HD1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HE2'})&&ismember('HE1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        
        % When chemical shift is given for one proton in a symmetric
        % aromatic ring, use it for the other one (usually true)
        pdb_chemsh{n}=bmrb_shifts(strcmp('HE1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                                                        pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
                                                    
    elseif ismember(pdb_atom_id{n},{'HH'})&&ismember(pdb_aa_typ(n),{'TYR'})
        
        % Ignore OH protons of tyrosines (usually rapid exchange)
        disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
        
    elseif ismember(pdb_atom_id{n},{'HD1'})&&ismember(pdb_aa_typ(n),{'HIS'})
        
        % Ignore ring NH protons of histidines (usually rapid exchange)
        disp(['WARNING: ring NH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
        
    elseif ismember(pdb_atom_id{n},{'NZ','HZ1','HZ2','HZ3'})&&ismember(pdb_aa_typ(n),{'LYS'})
        
        % Ignore side chain amino protons of lysines (usually rapid exchange)
        disp(['WARNING: side chain NH3+ group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
        
    elseif ismember(pdb_atom_id{n},{'NE','CZ','NH1','NH2','HE','HH11','HH12','HH21','HH22'})&&ismember(pdb_aa_typ(n),{'ARG'})
        
        % Ignore side chain amino protons of arginines (usually rapid exchange)
        disp(['WARNING: side chain nitrogen block of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
        
    end
    
end

% Estimate J-couplings
scalar_couplings=guess_j_pro(pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords);

% Estimate chemical shielding anisotropies
CSAs=guess_csa_pro(pdb_aa_num,pdb_atom_id,pdb_coords);

% Assign isotopes and labels
isotopes=cell(1,numel(pdb_atom_id));
for n=1:numel(pdb_atom_id) 
    switch pdb_atom_id{n}(1)
        case 'H'
            isotopes{n}='1H';
        case 'C'
            isotopes{n}='13C';
        case 'N'
            isotopes{n}='15N';
        otherwise
            error('unknown atom type.');
    end
end

% Process atom selection specification
if isnumeric(options.select)
    
    % Import atoms with user-specified numbers
    subset=false(size(pdb_atom_id)); subset(options.select)=true;

elseif strcmp(options.select,'backbone')
    
    % Import backbone up to CB and HB
    subset=ismember(pdb_atom_id,{'H','N','C','CA','HA','HA2','HA3',...
                                 'CB','HB','HB1','HB2','HB3'});
    
elseif strcmp(options.select,'backbone-minimal')
    
    % Import minimal backbone
    subset=ismember(pdb_atom_id,{'H','N','C','CA','HA','HA2','HA3'});
    
elseif strcmp(options.select,'backbone-hsqc')
    
    % Import backbone up to CB and HB and also ASN and GLN amide groups
    subset=ismember(pdb_atom_id,{'H','N','C','CA','HA','HA2','HA3','NE2',...
                                 'HE21','HE22','CD','CG','ND2','HD21',...
                                 'HD22','CB','HB','HB1','HB2','HB3'});
                             
elseif strcmp(options.select,'all')
    
    % Import everything
    subset=true(size(pdb_atom_id));

else 
    
    % Complain and bomb out
    error('incorrect subset selection specification.');
    
end

% Find missing chemical shifts
missing_shifts=find(cellfun(@isempty,pdb_chemsh))';
disp(' '); disp('############# SUMMARY OF MISSING ASSIGNMENTS ###############');

% Process unassigned chemical shifts
if strcmp(options.noshift,'keep')
    
    % Put unassigned spins between -1.0 and 0.0 ppm
    erzatz_shifts=linspace(-1,0,numel(missing_shifts));
    for n=1:numel(missing_shifts)
        disp(['Missing assignment of ' pdb_aa_typ{missing_shifts(n)} '(' num2str(pdb_aa_num(missing_shifts(n)))...
              ')-' pdb_atom_id{missing_shifts(n)} ': the spin is placed at ' num2str(erzatz_shifts(n)) ' ppm.']);
        pdb_chemsh{missing_shifts(n)}=erzatz_shifts(n);
    end
    
elseif strcmp(options.noshift,'delete')
    
    % Delete unassigned spins
    for n=missing_shifts
        disp(['Missing assignment of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ...
              ': the atom will not appear in the simulation.']);
    end
    subset=subset&(~cellfun(@isempty,pdb_chemsh));
    
else
    
    % Complain and bomb out
    error('incorrect value of options.noshift parameter.');
    
end

% Deuterate specified protons
if ~isempty(options.deuterate)

    % Get the gamma ratio
    gamma_ratio=spin('2H')/spin('1H');
    
    % Treat specified protons as deuterons
    for n=1:numel(isotopes)
       if strcmp(isotopes{n},'1H')&&...
          ismember(pdb_atom_id{n},options.deuterate)
      
            % Assign the isotope
            isotopes{n}='2H';
            
            % Scale the J-couplings
            for k=1:numel(isotopes)
                scalar_couplings{k,n}=gamma_ratio*scalar_couplings{k,n};
                scalar_couplings{n,k}=gamma_ratio*scalar_couplings{n,k};
            end
            
            % Inform the user
            disp(['Proton ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') ...
                  ' replaced by deuterium. Shifts kept, J-couplings scaled down.']);
              
            % Update the label
            pdb_atom_id{n}(1)='D';
              
       end
    end
   
end

% Apply the selection: system
sys.isotopes=isotopes(subset);
sys.labels=pdb_atom_id(subset)';

% Apply the selection: interactions
inter.zeeman.scalar=pdb_chemsh(subset)';
inter.zeeman.matrix=CSAs(subset)';
inter.coupling.scalar=scalar_couplings(subset,subset);
inter.coordinates=pdb_coords(subset);

% Apply the selection: auxiliaries
aux.pdb_aa_num=pdb_aa_num(subset)';
aux.pdb_aa_typ=pdb_aa_typ(subset)';

end

% Consistency enforcement
function grumble(pdb_file,bmrb_file,options)
if ~ischar(pdb_file)
    error('pdb_file must be a character string specifying a file name.');
end
if ~ischar(bmrb_file)
    error('bmrb_file must be a character string specifying a file name.');
end
if ~isfield(options,'select')
    error('options.select switch must be specfied.');
elseif (~isnumeric(options.select))&&(~ismember(options.select,{'all','backbone','backbone-hsqc','backbone-minimal','backbone-extended'}))
    error('invalid value for options.select, please refer to the manual.');
end
if ~isfield(options,'pdb_mol')
    error('options.pdb_mol switch must be specfied.');
elseif (~isnumeric(options.pdb_mol))
    error('invalid value for options.pdb_mol, please refer to the manual.');
end
if ~isfield(options,'noshift')
    error('options.noshift switch must be specfied.');
elseif (~isnumeric(options.noshift))&&(~ismember(options.noshift,{'keep','delete'}))
    error('invalid value for options.noshift, please refer to the manual.');
end
end

% We saw that we'd been given a law to live by, a moral law, they called
% it, which punished those who observed it -- for observing it. The more 
% you tried to live up to it, the more you suffered; the more you cheated
% it, the bigger reward you got.
%
% Ayn Rand, "Atlas Shrugged"

