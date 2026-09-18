% Parser for .magres files written by CASTEP and other codes in the
% CCP-NC magres v1.0 format. Reads the [atoms] and [magres] blocks
% and returns the geometry and the magnetic resonance tensors, keyed
% to the atoms in the order in which the atom records appear in the
% file. Syntax:
%
%                     props=c2spinach(file_name)
%
% Parameters:
%
%    file_name  - the name of the *.magres file, a
%                 character string
%
% Outputs:
%
%   props.filename         - log file name
%
%   props.symbols          - atomic symbols, 1 x natoms cell
%
%   props.std_geom         - atomic coordinates, natoms x 3, Angstrom
%
%   props.natoms           - number of atoms
%
%   props.cst              - chemical shielding tensors relative to
%                            the bare nucleus in vacuum, ppm, 1 x
%                            natoms cell, printed component order:
%                            xx xy xz yx yy yz zx zy zz
%
%   props.efg              - EFG tensors, a.u., 1 x natoms cell
%
%   props.k_couplings      - isotropic reduced spin-spin couplings,
%                            natoms x natoms, in the units used by
%                            gparse.m (magres K times mu_N^2/h, Hz);
%                            g2spinach.m converts them into J-coup-
%                            lings for the isotopes it is given
%
% Only the tensors that the file contains are returned; the caller
% should test for their presence with isfield. An atom for which a
% tensor is not printed gets an empty cell. Self-coupling records
% are skipped, so that the diagonal of props.k_couplings is zero
% as in gparse.m; when a file has both K(A,B) and K(B,A) records,
% their isotropic parts are averaged.
%
% Notes: CASTEP glues the atom label and the atom index into one
%        token when the index has three or more digits, for example
%        "ms O100"; such records are resolved against the atom list,
%        and a record that matches no atom or several atoms is an
%        error. Units records are checked against the standard ones
%        (Angstrom, ppm, au, 10^19.T^2.J^-1) and anything else is
%        an error.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=c2spinach.m>

function props=c2spinach(file_name)

% Check consistency
grumble(file_name);

% Read the file, drop the comments, and trim the lines
file_id=fopen(file_name,'r');
magres_log=textscan(file_id,'%s','delimiter','\n','whitespace','');
fclose(file_id); magres_log=strtrim(regexprep(magres_log{1},'#.*$',''));
props.filename=file_name;

% Locate the block tags and refuse nested or unbalanced blocks
tag_tokens=regexp(magres_log,'^[\[<](/?)([A-Za-z_]\w*)[\]>]$','tokens','once');
tag_lines=find(~cellfun(@isempty,tag_tokens)); open_line=0;
block_names={}; block_first=[]; block_last=[];
for n=tag_lines'
    if isempty(tag_tokens{n}{1})&&(open_line==0)
        open_line=n;
    elseif (open_line>0)&&strcmp(tag_tokens{n}{2},tag_tokens{open_line}{2})
        block_names{end+1}=tag_tokens{n}{2}; %#ok<AGROW>
        block_first(end+1)=open_line+1; block_last(end+1)=n-1; open_line=0; %#ok<AGROW>
    else
        error(['unexpected block tag in line ' num2str(n) ': ' magres_log{n}]);
    end
end
if open_line>0, error(['block [' tag_tokens{open_line}{2} '] is never closed.']); end
atoms_idx=find(strcmp(block_names,'atoms')); magres_idx=find(strcmp(block_names,'magres'));
if numel(atoms_idx)~=1, error('the file must contain exactly one [atoms] block.'); end
if numel(magres_idx)~=1, error('the file must contain exactly one [magres] block.'); end
atoms_block=magres_log(block_first(atoms_idx):block_last(atoms_idx));
magres_block=magres_log(block_first(magres_idx):block_last(magres_idx));

% Check the units of the records that are read below
units_tokens=regexp([atoms_block; magres_block],'^units\s+(\S+)\s+(\S+)$','tokens','once');
units_tokens=vertcat(units_tokens{:}); unit_tags={'atom','ms','efg','isc'};
unit_vals={'Angstrom','ppm','au','10^19.T^2.J^-1'};
for n=1:size(units_tokens,1)
    tag_idx=find(strcmp(unit_tags,units_tokens{n,1}),1);
    if ~isempty(tag_idx)&&~strcmp(units_tokens{n,2},unit_vals{tag_idx})
        error(['unsupported units for ' units_tokens{n,1} ' records: ' units_tokens{n,2}]);
    end
end

% Regular expressions for a number and for a 3x3 tensor
num_pat='([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)'; ten_pat=repmat(['\s+' num_pat],1,9);

% Atom records: species, label, index in label, and coordinates
atom_lines=find(~cellfun(@isempty,regexp(atoms_block,'^atom\s','once')));
atom_tokens=regexp(atoms_block(atom_lines),['^atom\s+(\S+)\s+(\S+)\s+(\d+)' repmat(['\s+' num_pat],1,3) '$'],'tokens','once');
bad_line=find(cellfun(@isempty,atom_tokens),1);
if ~isempty(bad_line), error(['malformed atom record: ' atoms_block{atom_lines(bad_line)}]); end
atom_tokens=vertcat(atom_tokens{:}); natoms=size(atom_tokens,1);
if natoms==0, error('the [atoms] block contains no atom records.'); end
props.symbols=atom_tokens(:,1)'; props.std_geom=str2double(atom_tokens(:,4:6));
props.natoms=natoms;

% Atom keys, delimited and glued, with canonical indices
atom_index=cellfun(@(x)num2str(str2double(x)),atom_tokens(:,3),'UniformOutput',false);
atom_keys=strcat(atom_tokens(:,2),{' '},atom_index)'; glued_keys=strcat(atom_tokens(:,2),atom_index)';
if numel(unique(atom_keys))~=natoms, error('duplicate atom label and index pairs in the [atoms] block.'); end

% Shielding and EFG records, matched to the atoms by label and index
cst=cell(1,natoms); efg=cell(1,natoms);
site_lines=find(~cellfun(@isempty,regexp(magres_block,'^(ms|efg)\s','once')));
site_tokens=regexp(magres_block(site_lines),['^(ms|efg)\s+(\S+(?:\s+\S+)?)' ten_pat '$'],'tokens','once');
for n=1:numel(site_lines)
    if isempty(site_tokens{n}), error(['malformed tensor record: ' magres_block{site_lines(n)}]); end
    atom_idx=key_atoms(strsplit(site_tokens{n}{2}),atom_keys,glued_keys,1);
    if size(atom_idx,1)~=1
        error(['tensor record does not match exactly one atom: ' magres_block{site_lines(n)}]);
    end
    tensor=reshape(str2double(site_tokens{n}(3:11)),[3 3])';
    if strcmp(site_tokens{n}{1},'ms')
        if ~isempty(cst{atom_idx}), error(['repeated ms record: ' magres_block{site_lines(n)}]); end
        cst{atom_idx}=tensor;
    else
        if ~isempty(efg{atom_idx}), error(['repeated efg record: ' magres_block{site_lines(n)}]); end
        efg{atom_idx}=tensor;
    end
end
if any(~cellfun(@isempty,cst)), props.cst=cst; end
if any(~cellfun(@isempty,efg)), props.efg=efg; end

% Reduced spin-spin coupling records, isotropic parts in gparse units
isc_lines=find(~cellfun(@isempty,regexp(magres_block,'^isc\s','once')));
isc_tokens=regexp(magres_block(isc_lines),['^isc\s+(\S+(?:\s+\S+){1,3})' ten_pat '$'],'tokens','once');
if ~isempty(isc_lines)
    k_sum=zeros(natoms,natoms); k_count=zeros(natoms,natoms);
    k_factor=1e19*(5.0507837461e-27)^2/6.62607015e-34;
    for n=1:numel(isc_lines)
        if isempty(isc_tokens{n}), error(['malformed isc record: ' magres_block{isc_lines(n)}]); end
        atom_pair=key_atoms(strsplit(isc_tokens{n}{1}),atom_keys,glued_keys,2);
        if size(atom_pair,1)~=1
            error(['isc record does not match exactly one atom pair: ' magres_block{isc_lines(n)}]);
        end
        if atom_pair(1)==atom_pair(2), continue; end
        k_iso=k_factor*sum(str2double(isc_tokens{n}([2 6 10])))/3;
        k_sum(atom_pair(1),atom_pair(2))=k_sum(atom_pair(1),atom_pair(2))+k_iso;
        k_sum(atom_pair(2),atom_pair(1))=k_sum(atom_pair(2),atom_pair(1))+k_iso;
        k_count(atom_pair(1),atom_pair(2))=k_count(atom_pair(1),atom_pair(2))+1;
        k_count(atom_pair(2),atom_pair(1))=k_count(atom_pair(2),atom_pair(1))+1;
    end
    props.k_couplings=k_sum./max(k_count,1);
end

end

% Resolves record label and index tokens into atom numbers, one row per viable split into delimited and glued keys
function atoms=key_atoms(id_tokens,atom_keys,glued_keys,natoms_needed)
splits={1,2,[1 1],[1 2],[2 1],[2 2]};
splits=splits(cellfun(@(x)(numel(x)==natoms_needed)&&(sum(x)==numel(id_tokens)),splits));
atoms=zeros(0,natoms_needed);
for n=1:numel(splits)
    candidate=zeros(1,natoms_needed); pos=1;
    for k=1:natoms_needed
        if splits{n}(k)==1
            match=find(strcmp(glued_keys,id_tokens{pos}));
        elseif all(isstrprop(id_tokens{pos+1},'digit'))
            match=find(strcmp(atom_keys,[id_tokens{pos} ' ' num2str(str2double(id_tokens{pos+1}))]));
        else
            match=[];
        end
        if numel(match)~=1, candidate=[]; break; end
        candidate(k)=match; pos=pos+splits{n}(k);
    end
    atoms=[atoms; candidate]; %#ok<AGROW>
end
end

% Consistency enforcement
function grumble(file_name)
if ~ischar(file_name)
    error('file_name must be a character string.');
end
end

% Malevolence lurks in men who avoid wine, games, the company
% of beautiful women, and conversations at dinner. Such people
% are either gravely ill, or they secretly hate those around.
%
% Mikhail Bulgakov, "Master and Margarita"

