% Reads BMRB file, extracts amino acid numbers, amino acid types,
% PDB atom identifiers and chemical shifts. Syntax:
%
%    [aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)
%
% Output parameters:
%
%       aa_num - amino acid numbers, vector
%
%       aa_typ - amino acid types, cell array of strings
%
%       pdb_id - atom identifiers using PDB nomenclature,
%                cell array of strings
%
%       chems  - chemical shifts in ppm
%
% Note: direct calls are discouraged, see the protein HOWTO docu-
%       ment for instructions on importing protein data.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=Read_bmrb.m>

function [aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)

% Check consistency
grumble(bmrb_file_name);

% Open the BMRB file
file_id=fopen(bmrb_file_name,'r');

% Get the outputs started
aa_num=[]; aa_typ={};
pdb_id={}; chemsh=[];

% Parse the BMRB file
while ~feof(file_id)
    
    % Get a line
    data_line=fgetl(file_id);
    
    % Replace tabs with spaces
    data_line=regexprep(data_line,'\t',' ');
    
    % If not empty, read
    if ~isempty(data_line)
        parsed_string=textscan(data_line,'%f %f %s %s %s %f %f %f','delimiter',' ','MultipleDelimsAsOne',1);
        if all(~cellfun(@isempty,parsed_string))
            aa_num(end+1)=parsed_string{2};    %#ok<AGROW>
            aa_typ{end+1}=parsed_string{3}{1}; %#ok<AGROW>
            pdb_id{end+1}=parsed_string{4}{1}; %#ok<AGROW>
            chemsh(end+1)=parsed_string{6};    %#ok<AGROW>
        end
    end
    
end

% Capitalize amino acid type specifications
for n=1:numel(aa_typ)
    aa_typ{n}=upper(aa_typ{n}); %#ok<AGROW>
end

% Make outputs column vectors
aa_num=aa_num'; aa_typ=aa_typ';
pdb_id=pdb_id'; chemsh=chemsh';

% Close the BMRB file
fclose(file_id);

% Check for empty returns
if isempty(chemsh)
    error('BMRB file could not be parsed - the data format is not canonical.');
end

end

% Consistency enforcement
function grumble(bmrb_file_name)
if ~ischar(bmrb_file_name)
    error('bmrb_file_name must be a character string.');
end
end

% If it's not true, it's well invented.
% 
% Dante

