% Merges multiple sys and inter structures into one. Useful for
% setting up chemical kinetics simulations where the molecules
% come from different DFT calculations. Syntax:
%
%         [sys,inter]=merge_inp(sys_parts,inter_parts)
%
% Parameters:
%
%    sys_parts   - a cell array of sys structures
%                  to be merged
%
%    inter_parts - a cell array of inter structures
%                  to be merged
%
% Outputs:
%
%    sys         - resulting sys structure
%
%    inter       - resulting inter structure
%
% Note: incomplete, an error is thrown for unhandled subfields,
%       extend as appropriate.
%
% ilya.kuprov@weizmann.ac.il
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=merge_inp.m>

function [sys,inter]=merge_inp(sys_parts,inter_parts)

% Check consistency
grumble(sys_parts,inter_parts); 

% Create structure stubs
sys.stub=1; inter.stub=1;

% Fields containing common scalars that are equal 
[sys,sys_parts]=merge_like_magnet(sys,sys_parts,'magnet');

% Fields containing row cell arrays
[sys,sys_parts]=merge_like_isotopes(sys,sys_parts,'isotopes');
[sys,sys_parts]=merge_like_isotopes(sys,sys_parts,'labels');

% Fields containing column cell arrays
[inter,inter_parts]=merge_like_coords(inter,inter_parts,'coordinates');

% Zeeman interaction specifications
zeeman_parts=strip(inter_parts,'zeeman'); zeeman.stub=1;
zeeman=merge_like_isotopes(zeeman,zeeman_parts,'matrix');
zeeman=merge_like_isotopes(zeeman,zeeman_parts,'scalar');
zeeman=merge_like_isotopes(zeeman,zeeman_parts,'eigs');
zeeman=merge_like_isotopes(zeeman,zeeman_parts,'euler');
zeeman=rmfield(zeeman,'stub'); inter.zeeman=zeeman;
inter_parts=cellfun(@(x)rmfield(x,'zeeman'),...
                    inter_parts,'UniformOutput',false);

% Coupling specifications
coupling_parts=strip(inter_parts,'coupling'); coupling.stub=1;
coupling=merge_like_couplings(coupling,coupling_parts,'matrix');
coupling=merge_like_couplings(coupling,coupling_parts,'scalar');
coupling=merge_like_couplings(coupling,coupling_parts,'eigs');
coupling=merge_like_couplings(coupling,coupling_parts,'euler');
coupling=rmfield(coupling,'stub'); inter.coupling=coupling;
inter_parts=cellfun(@(x)rmfield(x,'coupling'),...
                    inter_parts,'UniformOutput',false);

% Catch unhandled subfields
for n=1:numel(sys_parts)
    if ~isempty(fieldnames(sys_parts{n}))
        error('unhandled subfield in sys.');
    end
end
for n=1:numel(inter_parts)
    if ~isempty(fieldnames(inter_parts{n}))
        error('unhandled subfield in inter.');
    end
end

% Delete the stubs
sys=rmfield(sys,'stub'); inter=rmfield(inter,'stub');

end

% Consistency enforcement
function grumble(sys_parts,inter_parts)
if (~iscell(sys_parts))||(~iscell(inter_parts))
    error('both inputs must be cell arrays of data structures.');
end
end

% Strip a common field from an array of structures
function parts_array=strip(parts_array,field_name)
    parts_array=cellfun(@(x)x.(field_name),parts_array,'UniformOutput',false);
end

% Merge subsystem fields assuming they hold an optional common scalar
function [spec,spec_parts]=merge_like_magnet(spec,spec_parts,field_name)
if isfield(spec_parts{1},field_name)
    spec.(field_name)=spec_parts{1}.(field_name);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if (~isfield(spec,field_name))||...
               (spec_parts{n}.(field_name)~=spec.(field_name))
                error(['values of ' field_name ' are not present or not the same in all subsystems.']);
            else
                spec_parts{n}=rmfield(spec_parts{n},field_name);
            end
        else
            if isfield(spec,field_name)
                error(['values of ' field_name ' are not present or not the same in all subsystems.']);
            end
        end
    end
else
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            error(['values of ' field_name ' are not present or not the same in all subsystems.']);
        end
    end
end
end

% Merge subsystem fields assuming optional row cell arrays
function [spec,spec_parts]=merge_like_isotopes(spec,spec_parts,field_name)
if isfield(spec_parts{1},field_name)
    spec.(field_name)=spec_parts{1}.(field_name);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if ~isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            else
                spec.(field_name)=[spec.(field_name) spec_parts{n}.(field_name)];
                spec_parts{n}=rmfield(spec_parts{n},field_name);
            end
        else
            if isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            end
        end
    end
else
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            error([field_name ' is not present in all subsystems.']);
        end
    end
end
end

% Merge subsystem fields assuming optional column cell arrays
function [spec,spec_parts]=merge_like_coords(spec,spec_parts,field_name)
if isfield(spec_parts{1},field_name)
    spec.(field_name)=spec_parts{1}.(field_name);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if ~isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            else
                spec.(field_name)=[spec.(field_name); spec_parts{n}.(field_name)];
                spec_parts{n}=rmfield(spec_parts{n},field_name);
            end
        else
            if isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            end
        end
    end
else
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            error([field_name ' is not present in all subsystems.']);
        end
    end
end
end

% Merge subsystem fields assuming optional square cell arrays
function [spec,spec_parts]=merge_like_couplings(spec,spec_parts,field_name)
if isfield(spec_parts{1},field_name)
    spec.(field_name)=spec_parts{1}.(field_name);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if ~isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            else
                spec.(field_name)=blkdiag(spec.(field_name),spec_parts{n}.(field_name));
                spec_parts{n}=rmfield(spec_parts{n},field_name);
            end
        else
            if isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            end
        end
    end
else
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            error([field_name ' is not present in all subsystems.']);
        end
    end
end
end

% Dicebant ergo Pilato pontifices Judaeorum: "Noli scribere
% 'Rex Judaeorum' sed quia ipse dixit 'Rex sum Judaeorum'."
% Respondit Pilatus: "Quod scripsi, scripsi."
%
% Ioannes 19:21-22

