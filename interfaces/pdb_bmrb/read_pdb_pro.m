% Reads the a PDB file and returns amino acid numbers, amino acid types,
% PDB atom identifiers and Cartesian coordinates. Syntax:
%
%  [aa_num,aa_typ,pdb_id,coords]=read_pdb_pro(pdb_file_name,instance)
%
% Parameters:
%
%    pdb_file_name - a character string with the file name
%
%    mod_id        - the number of model that should be
%                    read from the pdb file
%
% Outputs:
%
%       aa_num - nspins x 1 vector giving the number of
%                the amino acid to which each spin belongs
%
%       aa_typ - nspins x 1 cell array of strings giving 
%                the PDB identifier of the amino acid to 
%                which each spin belongs (e.g. 'TYR')
%
%       pdb_id - nspins x 1 cell array of strings giving 
%                the PDB identifier of the protein atom 
%                type to which each spin belongs (e.g. 'HE2')
%
%       coords - nspins x 1 cell array of 3-vectors giving 
%                Cartesian coordinates of each spin in Angstrom
%
% i.ka.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=read_pdb_pro.m>

function [aa_num,aa_typ,pdb_id,coords]=read_pdb_pro(pdb_file_name,mod_id)

% Check consistency
grumble(pdb_file_name,mod_id);

% Open the PDB file
file_id=fopen(pdb_file_name,'r');

% Scroll to the selected structure
while ~feof(file_id)
    parsed_string=textscan(fgetl(file_id),'MODEL %d','delimiter',' ','MultipleDelimsAsOne',1);
    if parsed_string{1}==mod_id
        disp(['Reading model ' num2str(mod_id) ' from ' pdb_file_name '...' ]); break;
    end
end

% Get the outputs started
aa_num=[]; aa_typ={};
pdb_id={}; coords={};

% Parse the PDB file
while ~feof(file_id)
    data_line=fgetl(file_id);
    if (numel(data_line)>=6)&&strcmp('ENDMDL',data_line(1:6)), break; end
    if ~isempty(data_line)
        parsed_string=textscan(data_line,'ATOM %f %s %s %s %f %f %f %f %f %f','delimiter',' ','MultipleDelimsAsOne',1); 
        if all(~cellfun(@isempty,parsed_string))
            aa_num(end+1)=parsed_string{5};      %#ok<AGROW>
            aa_typ{end+1}=parsed_string{3}{1};   %#ok<AGROW>
            pdb_id{end+1}=parsed_string{2}{1};   %#ok<AGROW>
            coords{end+1}=[parsed_string{6:8}];  %#ok<AGROW>
        end
    end
end

% Capitalize amino acid type specifications
for n=1:numel(aa_typ)
    aa_typ{n}=upper(aa_typ{n}); %#ok<AGROW>
end

% Make outputs column vectors
aa_num=aa_num'; aa_typ=aa_typ';
pdb_id=pdb_id'; coords=coords';

% Close the PDB file
fclose(file_id);

end

% Consistency enforcement
function grumble(pdb_file_name,mod_id)
if ~ischar(pdb_file_name)
    error('pdb_file_name must be a character string.');
end
if (~isnumeric(mod_id))||(~isreal(mod_id))||...
   (mod_id<1)||(mod(mod_id,1)~=0)
    error('mod_id must be a positive real integer.');
end
end

% Heartless sterility, obliteration of all melody, all tonal charm, all
% music... This revelling in the destruction of all tonal essence, raging
% satanic fury in the orchestra, this demoniacal, lewd caterwauling,
% scandal-mongering, gun-toting music, with an orchestral accompaniment
% slapping you in the face... Hence, the secret fascination that makes it
% the darling of feeble-minded royalty... of the court monkeys covered with
% reptilian slime, and of the blase hysterical female court parasites who
% need this galvanic stimulation by massive instrumental treatment to throw
% their pleasure-weary frog-legs into violent convulsion... the diabolical
% din of this pig-headed man, stuffed with brass and sawdust, inflated, in
% an insanely destructive self-aggrandizement, by Mephistopheles' mephitic
% and most venomous hellish miasma, into Beelzebub's Court Composer and
% General Director of Hell's Music - Wagner!
%
% J.L. Klein, "Geschichte des Dramas", 1871.

