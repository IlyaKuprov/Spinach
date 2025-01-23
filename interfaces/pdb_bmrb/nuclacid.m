% Nucleic acid data import function. Parses PDB and chemical shift
% data, runs a J-coupling guess using guess_j_nuc.m function and 
% outputs sys and inter data structures that are required by create.m 
% gateway function in Spinach. Syntax:
%
%          [sys,inter]=nuclacid(pdb_file,shift_file,options)
%
% Parameters:
%
%           pdb_file - a character string containing the name 
%                      of the PDB file
%   
%         shift_file - a character string containing the name 
%                      of the chemical shift file, ASCII for-
%                      matted as [residue_number atom_id shift],
%                      see example.txt in examples/nmr_nucleic 
%      
%  options.deut_list - a cell array of strings, specifying which
%                      atoms should be assumed to be deuterated,
%                      for example {'ADE:H2pp'}. When an atom is
%                      deuterated, J-couplings are reduced appro-
%                      priately.
%  
%    options.noshift - 'keep' places unassigned atoms between -1 
%                       and 0 ppm, 'delete' removes them from the
%                       system
%
% 
% Returns:
%
%
%    sys.isotopes          - Nspins x 1 cell array of strings
%
%    sys.labels            - Nspins x 1 cell array of strings 
%                            containing standard IUPAC DNA/RNA 
%                            atom labels
%
%    inter.coordinates     - Nspins x 3 matrix, Angstrom.
%
%    inter.zeeman.scalar   - Nspins x 1 cell array of numbers, 
%                            ppm. Isotropic chemical shifts go
%                            here.
%
%    inter.coupling.scalar - Nspins x Nspins cell array of sca-
%                            lar couplings, all in Hz.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=Nuclacid.m>

function [sys,inter]=nuclacid(pdb_file,shift_file,options)

% Check consistency
grumble(pdb_file,shift_file,options);

% Parse the PDB file
[pdb_res_num,pdb_res_typ,pdb_atom_id,pdb_coords]=read_pdb_nuc(pdb_file);

% Open the shift file
file_id=fopen(shift_file,'r');

% Get the outputs started
bmrb_res_num=[]; bmrb_atom_id={}; bmrb_chemsh={};

% Parse the shift file
while ~feof(file_id)
    data_line=fgetl(file_id);
    if ~isempty(data_line)
        parsed_string=textscan(data_line,'%f %s %f','delimiter','\t','MultipleDelimsAsOne',1);
        if all(~cellfun(@isempty,parsed_string))
            bmrb_res_num(end+1)=parsed_string{1};    %#ok<AGROW>
            bmrb_atom_id{end+1}=parsed_string{2}{1}; %#ok<AGROW>
            bmrb_chemsh{end+1}=parsed_string{3};     %#ok<AGROW>
        end
    end
end

% Make outputs column vectors
bmrb_res_num=bmrb_res_num';
bmrb_atom_id=bmrb_atom_id';
bmrb_chemsh=bmrb_chemsh';

% Close the shift file
fclose(file_id);

% Remove oxygen, phosphorus and the exchangeable protons
kill_mask=ismember(pdb_atom_id,{'O5''','O4''','O4','O6','O2','O2''','O3''','O1P','O2P','P','H2''','H5T','HO''2'});
pdb_res_num(kill_mask)=[]; pdb_atom_id(kill_mask)=[]; 
pdb_res_typ(kill_mask)=[]; pdb_coords(kill_mask)=[];
disp('WARNING: oxygen atoms will not appear in the simulation.');

% Match chemical shifts
pdb_chemsh=cell(numel(pdb_atom_id),1);
disp(' '); disp('############# CHEMICAL SHIFT IMPORT SUMMARY ###############');
for n=1:numel(pdb_atom_id)
    
    % Pull the current amino acid from BMRB
    select_mask=(bmrb_res_num==pdb_res_num(n));
    bmrb_atoms=bmrb_atom_id(select_mask);
    bmrb_shifts=bmrb_chemsh(select_mask);
    
    % Match PDB and BMRB data
    if ismember(pdb_atom_id{n},bmrb_atoms)
        pdb_chemsh{n}=bmrb_shifts(strcmp(pdb_atom_id{n},bmrb_atoms)); pdb_chemsh{n}=pdb_chemsh{n}{1};
        disp(['Imported chemical shift for ' pad([pdb_res_typ{n} '(' num2str(pdb_res_num(n)) '):' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    end
    
end

% Replace primes with 'p' symbols for convenience
for n=1:numel(pdb_atom_id)
    pdb_atom_id{n}=strrep(pdb_atom_id{n},'''','p');
end

% Estimate J-couplings
scalar_couplings=guess_j_nuc(pdb_res_num,pdb_res_typ,pdb_atom_id,pdb_coords);

% CSA estimation goes here

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
        case 'P'
            isotopes{n}='31P';
        otherwise
            error('unknown atom type.');
    end
end

% Find missing chemical shifts
missing_shifts=find(cellfun(@isempty,pdb_chemsh))';
disp(' '); disp('############# SUMMARY OF MISSING ASSIGNMENTS ###############');

% Process unassigned chemical shifts
subset=true(size(pdb_atom_id));
if strcmp(options.noshift,'keep')
    
    % Put unassigned spins between -1.0 and 0.0 ppm
    erzatz_shifts=linspace(-1,0,numel(missing_shifts));
    for n=1:numel(missing_shifts)
        disp(['Missing assignment of ' pdb_res_typ{missing_shifts(n)} '(' num2str(pdb_res_num(missing_shifts(n)))...
              ')-' pdb_atom_id{missing_shifts(n)} ': the spin is placed at ' num2str(erzatz_shifts(n)) ' ppm.']);
        pdb_chemsh{missing_shifts(n)}=erzatz_shifts(n);
    end
    
elseif strcmp(options.noshift,'delete')
    
    % Delete unassigned spins
    for n=missing_shifts
        disp(['Missing assignment of ' pdb_res_typ{n} '(' num2str(pdb_res_num(n)) ')-' pdb_atom_id{n} ...
              ': the atom will not appear in the simulation.']);
    end
    subset=subset&(~cellfun(@isempty,pdb_chemsh));
    
else
    
    % Complain and bomb out
    error('incorrect value of options.noshift parameter.');
    
end
 
% Apply the selection
sys.isotopes=isotopes(subset);
inter.zeeman.scalar=pdb_chemsh(subset)';
inter.coupling.scalar=scalar_couplings(subset,subset);
inter.coordinates=pdb_coords(subset);

% Write detailed labels
pdb_atom_id=pdb_atom_id(subset);
pdb_res_num=pdb_res_num(subset);
pdb_res_typ=pdb_res_typ(subset);
sys.labels=cell(1,numel(pdb_atom_id));
for n=1:numel(pdb_atom_id)
    sys.labels{n}=[pdb_res_typ{n} '(' num2str(pdb_res_num(n)) '):' pdb_atom_id{n}];
end

% Deuterate the atoms from the list supplied
disp(' '); disp('############# SUMMARY OF ISOTOPE REASSIGNMENTS ###############');
for n=1:numel(sys.labels)
    if ismember([pdb_res_typ{n} ':' pdb_atom_id{n}],options.deut_list)
        if ~strcmp(sys.isotopes{n},'1H')
            error('the particle you are trying to deuterate is not a proton.');
        else
            sys.isotopes{n}='2H';
            disp(['Isotope type for ' sys.labels{n} ' changed to deuterium, J-coupling guess modified accordingly.']);
            inter.coupling.scalar(n,:)=inter.coupling.scalar(n,:)*(spin('2H')/spin('1H'));
            inter.coupling.scalar(:,n)=inter.coupling.scalar(:,n)*(spin('2H')/spin('1H'));
        end
    end
end

end

% Consistency enforcement
function grumble(pdb_file,shift_file,options)
if ~ischar(pdb_file)
    error('pdb_file must be a character string specifying a file name.');
end
if ~ischar(shift_file)
    error('shift_file must be a character string specifying a file name.');
end
if ~isfield(options,'noshift')
    error('options.noshift switch must be specfied.');
elseif ~ischar(options.noshift)
    error('options.noshift must be a character string.');
elseif ~ismember(options.noshift,{'keep','delete'})
    error('invalid value for options.noshift, please refer to the manual.');
end
if ~isfield(options,'deut_list')
    error('options.deut_list must be supplied.');
elseif ~iscell(options.deut_list)
    error('options.deut_list must be a cell array of character strings.');
end
end

% All money is a matter of belief.
%
% Adam Smith

