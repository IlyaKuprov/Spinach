% Reads the coordinates of all atoms from the user-specified PDB file
% and returns, for each atom, the residue number, the residue type, the
% PDB label and the Cartesian coordinates. Syntax:
%
%     [res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)
%
% Output parameters:
%
%     nuc_num  - nspins x 1 vector giving the number of the
%                nucleotide to which each spin belongs
% 
%     nuc_typ  - nspins x 1 cell array of strings giving the
%                PDB identifier of the nucleotide to which
%                each spin belongs (e.g. 'GUA')
% 
%     pdb_id   - nspins x 1 cell array of strings giving the
%                PDB identifier of the nucleic acid atom 
%                type to which each spin belongs (e.g. 'C1P')
% 
%     coords   - nspins x 1 cell array of 3-vectors giving 
%                Cartesian coordinates of each spin in Angstrom
%
% Note: All atoms in the file are read, make sure the PDB only contains
%       one model.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Read_pdb_nuc.m>

function [res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)

% Check consistency
grumble(pdb_file_name);

% Open the PDB file
file_id=fopen(pdb_file_name,'r');

% Get the outputs started
res_num=[]; res_typ={};
pdb_id={}; coords={};

% Parse the PDB file
while ~feof(file_id)
    data_line=fgetl(file_id);
    if ~isempty(data_line)
        parsed_string=textscan(data_line,'ATOM %f %s %s %f %f %f %f %f %f %s','delimiter',' ','MultipleDelimsAsOne',1); 
        if all(~cellfun(@isempty,parsed_string))
            res_num(end+1)=parsed_string{4};      %#ok<AGROW>
            res_typ{end+1}=parsed_string{3}{1};   %#ok<AGROW>
            pdb_id{end+1}=parsed_string{2}{1};    %#ok<AGROW>
            coords{end+1}=[parsed_string{5:7}];   %#ok<AGROW>
        end
    end
end

% Capitalise amino acid type specifications
for n=1:numel(res_typ)
    res_typ{n}=upper(res_typ{n}); %#ok<AGROW>
end

% Make outputs column vectors
res_num=res_num'; res_typ=res_typ';
pdb_id=pdb_id'; coords=coords';

% Close the PDB file
fclose(file_id);

end

% Consistency enforcement
function grumble(pdb_file_name)
if ~ischar(pdb_file_name)
    error('pdb_file_name must be a character string.');
end
end

% I think that it might well be premature to make an award of
% a Prize to Watson and Crick, because of existing uncertainty
% about the detailed structure of nucleic acid. I myself feel
% that it is likely that the general nature of the Watson-Crick
% structure is correct, but that there is doubt about details.
% With respect to Wilkins, I may say that I recognize his virtu-
% osity in having grown better fibers of DNA than any that had
% been grown before and in having obtained [better] x-ray photo-
% graphs than were available before, but I doubt that this work
% represents a sufficient contribution to chemistry to permit 
% him to be included among recipients of a Nobel Prize.
%
% Linus Pauling, in his letter to the Nobel Prize
% Committee, 15 March 1960.

