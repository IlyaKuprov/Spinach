% Tests tagged NMR-STAR chemical shift import and legacy tables. Syntax:
%
%                       result=test_bmrb_parser()
%
% Outputs:
%
%     result - checks for tag order, token boundaries, nulls,
%              author numbering, multiple loops, and malformed rows
%
% talos@spindynamics.org

function result=test_bmrb_parser()

% Describe the import contract
result=new_test_result('interfaces/bmrb_parser','BMRB shift parser',...
                       'tagged assignments retain their values independently of layout.');

% Create disposable files for production calls
work_dir=tempname; mkdir(work_dir);
cleanup=onCleanup(@()rmdir(work_dir,'s'));
file_name=fullfile(work_dir,'fixture.str');

% Exercise reordered columns, comments, quotes, and wrapped records
header=["data_test";"save_shifts";"loop_";...
        "_Atom_chem_shift.Val";"_Atom_chem_shift.Atom_ID";...
        "_Atom_chem_shift.Comp_ID";"_Atom_chem_shift.Comp_index_ID"];
lines=[header;"4.2 'HA2' gly 1 # first assignment";...
       "55.1";"CA GLY 1";"stop_";"save_"];
writelines(lines,file_name);
result=check_import(result,file_name,[1;1],{'GLY';'GLY'},...
                    {'HA2';'CA'},[4.2;55.1],'reordered and wrapped');

% Preserve apostrophes and hashes in values rather than treating them as syntax
lines=[header;"4.2 H5' G 1";"55.1 'C#1' G 1";"stop_";"save_"];
writelines(lines,file_name);
result=check_import(result,file_name,[1;1],{'G';'G'},...
                    {'H5''';'C#1'},[4.2;55.1],'punctuation');

% Ignore keywords in quoted and multiline metadata while collecting every list
lines=["data_test";"_Entry.Details";";";"loop_ _Atom_chem_shift.Val";...
       "1 2 3 4 unrelated numeric text";";";"_Entry.Title 'loop_'";...
       header(2:end);"4.2 HA2 GLY 1";"stop_";"save_";...
       header(2:end);"55.1 CA GLY 1";"stop_";"save_"];
writelines(lines,file_name);
result=check_import(result,file_name,[1;1],{'GLY';'GLY'},...
                    {'HA2';'CA'},[4.2;55.1],'metadata and lists');

% Exclude the opening text-field newline while retaining its data newlines
lines=[header;"4.2";";";"HA2";";";"GLY 1";"stop_";"save_"];
writelines(lines,file_name);
result=check_import(result,file_name,1,{'GLY'},{sprintf('HA2\n')},4.2,'text field value');

% Allow null optional fields and omit null assignments without synthesising values
lines=[header;". HA2 GLY 1";"? CA GLY 1";"4.2 HA2 GLY 1";"stop_";"save_"];
writelines(lines,file_name);
result=check_import(result,file_name,1,{'GLY'},{'HA2'},4.2,'null shifts');

% Match author residue numbers and fall back independently for each row
lines=[header;"_Atom_chem_shift.Auth_seq_ID";"_Atom_chem_shift.Seq_ID";...
       "4.2 HA2 GLY 1 -2 7";"55.1 CA GLY 1 . 8";...
       "8.1 H GLY 3 ? .";"stop_";"save_"];
writelines(lines,file_name);
result=check_import(result,file_name,[-2;8;3],{'GLY';'GLY';'GLY'},...
                    {'HA2';'CA';'H'},[4.2;55.1;8.1],'author fallback');

% Omit null residue assignments and preserve quoted structural words as labels
lines=[header;"4.2 'loop_' GLY 1";"4.3 HA2 GLY .";"stop_";"save_"];
writelines(lines,file_name);
result=check_import(result,file_name,1,{'GLY'},{'loop_'},4.2,'null residues');

% Read NMR-STAR 2.1 assignments by their tags rather than fixed positions
lines=["loop_";"_Atom_name";"_Chem_shift_value";...
       "_Residue_label";"_Residue_seq_code";"CA 55.1 gly 1";...
       "HA2 4.2 gly 1";"stop_"];
writelines(lines,file_name);
result=check_import(result,file_name,[1;1],{'GLY';'GLY'},...
                    {'CA';'HA2'},[55.1;4.2],'STAR 2.1');

% Ignore the similarly named chemical shift reference loop in STAR 2.1
lines=["data_test";"loop_";"_Mol_common_name";"_Chem_shift_value";...
       "TSP 0.0";"stop_";lines];
writelines(lines,file_name);
result=check_import(result,file_name,[1;1],{'GLY';'GLY'},...
                    {'CA';'HA2'},[55.1;4.2],'reference loop');

% Accept omitted stop markers at end of file and CRLF line endings
text=strjoin([header(1);header(3:end);"4.2 HA2 GLY 1"],sprintf('\r\n'));
file_id=fopen(file_name,'w'); fprintf(file_id,'%s',text); fclose(file_id);
result=check_import(result,file_name,1,{'GLY'},{'HA2'},4.2,'EOF and CRLF');

% Keep the historical headerless eight-column interpretation
writelines(["1 2 GLY CA C 55.1 . 1";"2 2 GLY HA2 H 4.2 0.02 1"],file_name);
result=check_import(result,file_name,2,{'GLY'},{'HA2'},4.2,'legacy layout');

% Return empty columns for valid files without shift assignments
writelines(["data_test";"save_metadata";"_Entry.Title 'no shifts'";"save_"],file_name);
result=check_import(result,file_name,zeros(0,1),cell(0,1),...
                    cell(0,1),zeros(0,1),'no assignments');

% Reject incomplete rows and invalid numerical assignments
invalids={[header;"4.2 HA2 GLY";"stop_";"save_"],...
          [header;"bad HA2 GLY 1";"stop_";"save_"],...
          [header;"4.2 HA2 GLY 1.5";"stop_";"save_"],...
          [header([1:3 5:end]);"HA2 GLY 1";"stop_";"save_"],...
          [header;header(4);"4.2 HA2 GLY 1 4.2";"stop_";"save_"],...
          [header;"""4.2 HA2 GLY 1";"stop_";"save_"],...
          [header;"4.2+1i HA2 GLY 1";"stop_";"save_"],...
          [header;"4.2 HA2 GLY 1+2i";"stop_";"save_"],...
          ["data_test";"_Entry.Details";";unterminated";header(2:end)]};
for n=1:numel(invalids)

    % Require a production error rather than truncated or fabricated data
    writelines(invalids{n},file_name); caught=false;
    try
        read_bmrb(file_name);
    catch
        caught=true;
    end
    result=test_true(result,['invalid table ' num2str(n)],caught,...
                     'malformed chemical shift assignments must be rejected');
end

end

% Compare every output of the production reader against the fixture
function result=check_import(result,file_name,nums,types,atoms,shifts,label)

% Read the fixture without duplicating parser logic
[aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(file_name);
result=test_true(result,label,isequal(aa_num,nums)&&isequal(aa_typ,types)&&...
                 isequal(pdb_id,atoms)&&isequal(chemsh,shifts),...
                 'all four outputs must match exactly, including their column shapes');

end


