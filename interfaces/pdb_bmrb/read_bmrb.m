% Reads BMRB file, extracts amino acid numbers, amino acid types,
% PDB atom identifiers, and chemical shifts. Syntax:
%
%    [aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)
%
% Parameters:
%
%       bmrb_file_name - a character string with the file name
%
% Outputs:
%
%       aa_num - amino acid numbers, vector
%
%       aa_typ - amino acid types, cell array of strings
%
%       pdb_id - atom identifiers using PDB nomenclature,
%                cell array of strings
%
%       chemsh - chemical shifts in ppm
%
% Supports tagged NMR-STAR 2.1 and 3.x shift loops, arbitrary column
% order, quoted and multiline values, comments, and wrapped rows.
% Numeric author residue numbers take precedence over Seq_ID and
% Comp_index_ID, matching PDB numbering where available. Missing
% author numbers fall back row by row; non-numeric author identifiers
% cannot be represented by this numeric interface. Null assignments
% are omitted. All shift lists and entities are returned in file order;
% this interface does not select between them. Valid STAR files without
% shift assignments return empty columns. Legacy eight-column tables
% remain supported. Atom labels are not renamed.
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

% Read and normalise line endings for STAR text fields
text=fileread(bmrb_file_name);
text=regexprep(text,'\r\n?','\n');

% Tokenise text fields, quoted values, comments, and bare values
pattern=['(?m)^;[\s\S]*?^;(?=\s|$)|'...
         '''(?:[^''\n]|''(?!\s|$))*''(?=\s|$)|'...
         '"(?:[^"\n]|"(?!\s|$))*"(?=\s|$)|'...
         '#[^\n]*|[^\s]+'];
tokens=regexp(text,pattern,'match');
tokens(startsWith(tokens,'#'))=[];
quoted=startsWith(tokens,{'''','"'})|...
       (startsWith(tokens,';')&contains(tokens,newline));

% Retain the established headerless eight-column import
if ~any(~quoted&(startsWith(lower(tokens),'data_')|...
                 strcmpi(tokens,'loop_')))

    % Open the BMRB file
    file_id=fopen(bmrb_file_name,'r');
    cleanup=onCleanup(@()fclose(file_id));

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

    % Finish the legacy outputs without changing their interpretation
    aa_typ=upper(aa_typ);
    aa_num=aa_num(:); aa_typ=aa_typ(:);
    pdb_id=pdb_id(:); chemsh=chemsh(:);
    if isempty(chemsh)
        error('BMRB file contains no recognised chemical shift table.');
    end
    return

end

% Distinguish structural keywords from quoted data values
lower_tokens=lower(tokens);
control=~quoted&(startsWith(lower_tokens,{'_','save_','data_'})|...
                  ismember(lower_tokens,{'loop_','stop_','global_'}));
for n=find(quoted)
    tokens{n}=tokens{n}(2:end-1);
end

% Collect the four documented assignment columns from every shift loop
aa_num=zeros(0,1); aa_typ=cell(0,1);
pdb_id=cell(0,1); chemsh=zeros(0,1);
for n=find(~quoted&strcmpi(lower_tokens,'loop_'))

    % Locate the tag list and the following loop values
    first=n+1; last=first;
    while last<=numel(tokens)&&~quoted(last)&&startsWith(tokens{last},'_')
        last=last+1;
    end
    tags=lower_tokens(first:last-1);
    finish=find(control(last:end),1);
    if isempty(finish)
        finish=numel(tokens)+1;
    else
        finish=last+finish-1;
    end

    % Resolve modern and legacy chemical shift tag names
    if any(startsWith(tags,'_atom_chem_shift.'))
        num_tags={'_atom_chem_shift.auth_seq_id',...
                  '_atom_chem_shift.seq_id','_atom_chem_shift.comp_index_id'};
        type_tag='_atom_chem_shift.comp_id';
        atom_tag='_atom_chem_shift.atom_id'; shift_tag='_atom_chem_shift.val';
    elseif any(strcmp(tags,'_atom_shift_assign_id'))||...
           all(ismember({'_residue_seq_code','_chem_shift_value'},tags))
        num_tags={'_residue_author_seq_code','_residue_seq_code'};
        type_tag='_residue_label'; atom_tag='_atom_name';
        shift_tag='_chem_shift_value';
    else
        continue
    end

    % Check the assignment table before interpreting its values
    num_cols=find(ismember(tags,num_tags));
    type_col=find(strcmp(tags,type_tag));
    atom_col=find(strcmp(tags,atom_tag));
    shift_col=find(strcmp(tags,shift_tag));
    if isempty(num_cols)||numel(type_col)~=1||...
       numel(atom_col)~=1||numel(shift_col)~=1||numel(unique(tags))~=numel(tags)
        error('Chemical shift loop has missing or duplicate assignment tags.');
    end
    if mod(finish-last,numel(tags))~=0
        error('Chemical shift loop has an incomplete row.');
    end
    table=reshape(tokens(last:finish-1),numel(tags),[])';

    % Prefer numeric author numbering with per-row sequence fallback
    residue=nan(size(table,1),1);
    for k=1:numel(num_tags)
        col=find(strcmp(tags,num_tags{k}),1);
        if ~isempty(col)
            values=str2double(table(:,col));
            mask=~isfinite(residue)&isfinite(values);
            residue(mask)=values(mask);
        end
    end

    % Omit null assignments without inventing numerical shifts
    shifts=str2double(table(:,shift_col));
    nulls={'.','?'};
    missing=all(ismember(table(:,num_cols),nulls),2)|...
            ismember(table(:,shift_col),nulls)|...
            ismember(table(:,type_col),nulls)|ismember(table(:,atom_col),nulls);
    if any(~missing&(~isfinite(shifts)|~isfinite(residue)|residue~=fix(residue)))
        error('Chemical shift assignment has an invalid shift or residue number.');
    end
    keep=~missing;
    aa_num=[aa_num;residue(keep)]; %#ok<AGROW>
    aa_typ=[aa_typ;upper(table(keep,type_col))]; %#ok<AGROW>
    pdb_id=[pdb_id;table(keep,atom_col)]; %#ok<AGROW>
    chemsh=[chemsh;shifts(keep)]; %#ok<AGROW>

end

end

% Consistency enforcement
function grumble(bmrb_file_name)
if ~ischar(bmrb_file_name)||~isrow(bmrb_file_name)||isempty(bmrb_file_name)
    error('bmrb_file_name must be a non-empty character row vector.');
end
end

% If it's not true, it's well invented.
%
% Dante


