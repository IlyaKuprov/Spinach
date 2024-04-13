% Reads the "new format" section of CASTEP .magres files. Syntax:
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
%   props.std_geom         - atomic coordinates (Angstrom)
%
%   props.symbols          - atomic symbols
%
%   props.natoms           - number of atoms
%
%   props.filename         - log file name
%
%   props.cst              - chemical shielding tensors relative
%                            to the bare nucleus in vacuum, ppm
%
%   props.efg              - EFG tensors, a.u.^-3
%
% Notes: CASTEP has a printing bug with over 100 nuclei. A kind of
%        workaround was implemented, but keep an eye on it.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=c2spinach.m>

function props=c2spinach(file_name)

% Check consistency
grumble(file_name);

% Read the file
file_id=fopen(file_name,'r');
castep_log=textscan(file_id,'%s','delimiter','\n');
fclose(file_id); castep_log=castep_log{1};
props.filename=file_name;

% Deblank all lines
for n=1:numel(castep_log)
    castep_log(n)=deblank(castep_log(n)); 
end

% Cartesian coordinates
props.symbols={}; props.std_geom=zeros(0,3);
for n=1:numel(castep_log)
    if (numel(castep_log{n})>3)&&strcmp(castep_log{n}(1:4),'atom')
        atom_spec=textscan(castep_log{n},'atom %s %s %f %f %f %f',...
                           'Delimiter',' ','MultipleDelimsAsOne',1);
        props.symbols{end+1}=atom_spec{1}{1};
        props.std_geom(end+1,:)=[atom_spec{4:6}];
    end
end

% Shielding tensors
props.cst={};
for n=1:numel(castep_log)
    if (numel(castep_log{n})>1)&&strcmp(castep_log{n}(1:2),'ms')
        cst_spec=textscan(castep_log{n},'ms %s %f %f %f %f %f %f %f %f %f %f',...
                          'Delimiter',' ','MultipleDelimsAsOne',1);
        if isempty(cst_spec{end}) % CASTEP printing bug
            cst_spec={cst_spec{1}{1}(1) str2double(cst_spec{1}{1}(2:end)) cst_spec{2:(end-1)}};
        end
        props.cst{end+1}=reshape([cst_spec{3:11}],[3 3]);
    end
end

% Electric field gradients
props.efg={};
for n=1:numel(castep_log)
    if (numel(castep_log{n})>2)&&strcmp(castep_log{n}(1:3),'efg')
        efg_spec=textscan(castep_log{n},'efg %s %f %f %f %f %f %f %f %f %f %f',...
                          'Delimiter',' ','MultipleDelimsAsOne',1);
        if isempty(efg_spec{end}) % CASTEP printing bug
            efg_spec={efg_spec{1}{1}(1) str2double(efg_spec{1}{1}(2:end)) efg_spec{2:(end-1)}};
        end
        props.efg{end+1}=reshape([efg_spec{3:11}],[3 3]);
    end
end

% Atom count
props.natoms=size(props.std_geom,1);
        
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

