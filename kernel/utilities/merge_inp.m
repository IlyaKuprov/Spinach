% Merges multiple sys and inter structures into one. Useful for
% setting up chemical kinetics simulations where the molecules
% come from different DFT calculations. Syntax:
%
%         [sys,inter]=merge_inp(sys_parts,inter_parts)
%
% Parameters:
%
%    sys_parts   - a row cell array of sys structures
%                  to be merged
%
%    inter_parts - a row cell array of inter structures
%                  to be merged
%
% Outputs:
%
%    sys         - resulting sys structure
%
%    inter       - resulting inter structure
%
% Note: extensive fields are concatenated with spin and chemical
%       subsystem indices offset as appropriate; non-extensive
%       fields (magnet, temperature, relaxation settings, etc.)
%       must be identical in all subsystems that supply them, and
%       any difference is treated as an error. Nested groups
%       (zeeman, coupling, giant, suscept, chem) and every field
%       must be present in all subsystems or in none. An error is
%       thrown for unhandled subfields. Coordinates and suscepti-
%       bility centres from all subsystems are assumed to refer
%       to one common frame of reference; spin index lists are
%       returned as row vectors.
%
% ilya.kuprov@weizmann.ac.il
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=merge_inp.m>

function [sys,inter]=merge_inp(sys_parts,inter_parts)

% Check consistency
grumble(sys_parts,inter_parts);

% Count the spins in each subsystem
part_sizes=cellfun(@(x)numel(x.isotopes),sys_parts);

% Create structure stubs
sys.stub=1; inter.stub=1;

% Fields containing common values that must be equal
sys_common={'magnet','output','scratch','enable',...
            'disable','tols','parallel','parprops'};
for n=1:numel(sys_common)
    [sys,sys_parts]=merge_like_magnet(sys,sys_parts,sys_common{n});
end

% Fields containing row cell arrays
[sys,sys_parts]=merge_like_isotopes(sys,sys_parts,'isotopes');
[sys,sys_parts]=merge_like_isotopes(sys,sys_parts,'labels');

% Fields containing common values that must be equal
inter_common={'relaxation','rlx_keep','rlx_dfs','equilibrium',...
              'temperature','damp_rate','srfk_tau_c',...
              'nz_shift','nz_onshell',...
              'weiz_r1e','weiz_r2e','weiz_r1n','weiz_r2n',...
              'nott_r1e','nott_r2e','nott_r1n','nott_r2n','pbc'};
for n=1:numel(inter_common)
    [inter,inter_parts]=merge_like_magnet(inter,inter_parts,inter_common{n});
end

% Fields that are per chemical species when a species split is present
if group_check(inter_parts,'chem')&&...
   all(cellfun(@(x)isfield(x.chem,'parts'),inter_parts))
    [inter,inter_parts]=merge_like_isotopes(inter,inter_parts,'tau_c');
    [inter,inter_parts]=merge_like_isotopes(inter,inter_parts,'order_matrix');
else
    [inter,inter_parts]=merge_like_magnet(inter,inter_parts,'tau_c');
    [inter,inter_parts]=merge_like_magnet(inter,inter_parts,'order_matrix');
end

% Fields containing column cell arrays
[inter,inter_parts]=merge_like_coords(inter,inter_parts,'coordinates');

% Fields containing per-spin relaxation rates
inter_rates={'r1_rates','r2_rates','lind_r1_rates','lind_r2_rates'};
for n=1:numel(inter_rates)
    [inter,inter_parts]=merge_like_rates(inter,inter_parts,inter_rates{n});
end

% Fields containing square arrays over spin pairs
inter_square={'srfk_mdepth','weiz_r1d','weiz_r2d'};
for n=1:numel(inter_square)
    [inter,inter_parts]=merge_like_couplings(inter,inter_parts,inter_square{n});
end

% Fields containing spin index vectors
[inter,inter_parts]=merge_like_sources(inter,inter_parts,'srsk_sources',part_sizes);

% Fields containing cell arrays of spin index vectors
[inter,inter_parts]=merge_like_parts(inter,inter_parts,'ignore',part_sizes);

% Zeeman interaction specifications
if group_check(inter_parts,'zeeman')
    zeeman_parts=strip(inter_parts,'zeeman'); zeeman.stub=1;
    [zeeman,zeeman_parts]=merge_like_isotopes(zeeman,zeeman_parts,'matrix');
    [zeeman,zeeman_parts]=merge_like_isotopes(zeeman,zeeman_parts,'scalar');
    [zeeman,zeeman_parts]=merge_like_isotopes(zeeman,zeeman_parts,'eigs');
    [zeeman,zeeman_parts]=merge_like_isotopes(zeeman,zeeman_parts,'euler');
    group_gate(zeeman_parts,'zeeman');
    zeeman=rmfield(zeeman,'stub'); inter.zeeman=zeeman;
    inter_parts=cellfun(@(x)rmfield(x,'zeeman'),...
                        inter_parts,'UniformOutput',false);
end

% Coupling specifications
if group_check(inter_parts,'coupling')
    coupling_parts=strip(inter_parts,'coupling'); coupling.stub=1;
    [coupling,coupling_parts]=merge_like_couplings(coupling,coupling_parts,'matrix');
    [coupling,coupling_parts]=merge_like_couplings(coupling,coupling_parts,'scalar');
    [coupling,coupling_parts]=merge_like_couplings(coupling,coupling_parts,'eigs');
    [coupling,coupling_parts]=merge_like_couplings(coupling,coupling_parts,'euler');
    group_gate(coupling_parts,'coupling');
    coupling=rmfield(coupling,'stub'); inter.coupling=coupling;
    inter_parts=cellfun(@(x)rmfield(x,'coupling'),...
                        inter_parts,'UniformOutput',false);
end

% Giant spin interaction specifications
if group_check(inter_parts,'giant')
    giant_parts=strip(inter_parts,'giant'); giant.stub=1;
    [giant,giant_parts]=merge_like_isotopes(giant,giant_parts,'coeff');
    [giant,giant_parts]=merge_like_isotopes(giant,giant_parts,'euler');
    group_gate(giant_parts,'giant');
    giant=rmfield(giant,'stub'); inter.giant=giant;
    inter_parts=cellfun(@(x)rmfield(x,'giant'),...
                        inter_parts,'UniformOutput',false);
end

% Magnetic susceptibility centre specifications
if group_check(inter_parts,'suscept')
    suscept_parts=strip(inter_parts,'suscept'); suscept.stub=1;
    [suscept,suscept_parts]=merge_like_isotopes(suscept,suscept_parts,'chi');
    [suscept,suscept_parts]=merge_like_isotopes(suscept,suscept_parts,'xyz');
    group_gate(suscept_parts,'suscept');
    suscept=rmfield(suscept,'stub'); inter.suscept=suscept;
    inter_parts=cellfun(@(x)rmfield(x,'suscept'),...
                        inter_parts,'UniformOutput',false);
end

% Chemical process specifications
if group_check(inter_parts,'chem')
    chem_parts=strip(inter_parts,'chem'); chem.stub=1;
    [chem,chem_parts]=merge_like_parts(chem,chem_parts,'parts',part_sizes);
    [chem,chem_parts]=merge_like_couplings(chem,chem_parts,'rates');
    [chem,chem_parts]=merge_like_isotopes(chem,chem_parts,'concs');
    [chem,chem_parts]=merge_like_couplings(chem,chem_parts,'flux_rate');
    [chem,chem_parts]=merge_like_magnet(chem,chem_parts,'flux_type');
    [chem,chem_parts]=merge_like_magnet(chem,chem_parts,'rp_theory');
    [chem,chem_parts]=merge_like_sources(chem,chem_parts,'rp_electrons',part_sizes);
    [chem,chem_parts]=merge_like_magnet(chem,chem_parts,'rp_rates');
    group_gate(chem_parts,'chem');
    chem=rmfield(chem,'stub'); inter.chem=chem;
    inter_parts=cellfun(@(x)rmfield(x,'chem'),...
                        inter_parts,'UniformOutput',false);
end

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

% Check that a nested group is present in all subsystems or none
function group_present=group_check(parts_array,group_name)
present=cellfun(@(x)isfield(x,group_name),parts_array);
if any(present)&&(~all(present))
    error([group_name ' specification must be present in all subsystems or none.']);
end
group_present=all(present);
end

% Catch unhandled subfields inside a nested group
function group_gate(group_parts,group_name)
for n=1:numel(group_parts)
    if ~isempty(fieldnames(group_parts{n}))
        error(['unhandled subfield in inter.' group_name '.']);
    end
end
end

% Strip a common field from an array of structures
function parts_array=strip(parts_array,field_name)
    parts_array=cellfun(@(x)x.(field_name),parts_array,'UniformOutput',false);
end

% Merge subsystem fields assuming they hold an optional common value
function [spec,spec_parts]=merge_like_magnet(spec,spec_parts,field_name)
if isfield(spec_parts{1},field_name)
    spec.(field_name)=spec_parts{1}.(field_name);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if (~isfield(spec,field_name))||...
               (~isequal(spec_parts{n}.(field_name),spec.(field_name)))
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

% Merge subsystem fields assuming optional per-spin rate arrays
function [spec,spec_parts]=merge_like_rates(spec,spec_parts,field_name)
if isfield(spec_parts{1},field_name)
    spec.(field_name)=spec_parts{1}.(field_name)(:);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if ~isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            else
                spec.(field_name)=[spec.(field_name); spec_parts{n}.(field_name)(:)];
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

% Merge subsystem fields assuming optional square arrays
function [spec,spec_parts]=merge_like_couplings(spec,spec_parts,field_name)
if isfield(spec_parts{1},field_name)
    spec.(field_name)=spec_parts{1}.(field_name);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if ~isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            elseif iscell(spec.(field_name))~=iscell(spec_parts{n}.(field_name))
                error([field_name ' must be of the same type in all subsystems.']);
            elseif iscell(spec.(field_name))
                old_block=spec.(field_name); new_block=spec_parts{n}.(field_name);
                merged=cell(size(old_block,1)+size(new_block,1));
                merged(1:size(old_block,1),1:size(old_block,1))=old_block;
                merged((size(old_block,1)+1):end,(size(old_block,1)+1):end)=new_block;
                spec.(field_name)=merged;
                spec_parts{n}=rmfield(spec_parts{n},field_name);
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

% Merge subsystem fields assuming optional spin index vectors
function [spec,spec_parts]=merge_like_sources(spec,spec_parts,field_name,part_sizes)
offsets=cumsum([0 part_sizes(1:end-1)]);
if isfield(spec_parts{1},field_name)
    spec.(field_name)=reshape(spec_parts{1}.(field_name),1,[]);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if ~isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            else
                spec.(field_name)=[spec.(field_name) ...
                                   reshape(spec_parts{n}.(field_name),1,[])+offsets(n)];
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

% Merge subsystem fields assuming optional cell arrays of spin index vectors
function [spec,spec_parts]=merge_like_parts(spec,spec_parts,field_name,part_sizes)
offsets=cumsum([0 part_sizes(1:end-1)]);
if isfield(spec_parts{1},field_name)
    spec.(field_name)=reshape(spec_parts{1}.(field_name),1,[]);
    spec_parts{1}=rmfield(spec_parts{1},field_name);
    for n=2:numel(spec_parts)
        if isfield(spec_parts{n},field_name)
            if ~isfield(spec,field_name)
                error([field_name ' is not present in all subsystems.']);
            else
                shifted=cellfun(@(x)x+offsets(n),...
                                reshape(spec_parts{n}.(field_name),1,[]),...
                                'UniformOutput',false);
                spec.(field_name)=[spec.(field_name) shifted];
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

% Consistency enforcement
function grumble(sys_parts,inter_parts)
if (~iscell(sys_parts))||(~iscell(inter_parts))
    error('both inputs must be cell arrays of data structures.');
end
if numel(sys_parts)~=numel(inter_parts)
    error('sys_parts and inter_parts must have the same number of elements.');
end
if isempty(sys_parts)
    error('at least one subsystem must be supplied.');
end
if (~isrow(sys_parts))||(~isrow(inter_parts))
    error('sys_parts and inter_parts must be row cell arrays.');
end
if any(cellfun(@(x)(~isstruct(x)),sys_parts))||...
   any(cellfun(@(x)(~isstruct(x)),inter_parts))
    error('all elements of sys_parts and inter_parts must be structures.');
end
if any(cellfun(@(x)(~isfield(x,'isotopes')),sys_parts))
    error('every sys structure must contain an isotopes subfield.');
end
end

% Dicebant ergo Pilato pontifices Judaeorum: "Noli scribere
% 'Rex Judaeorum' sed quia ipse dixit 'Rex sum Judaeorum'."
% Respondit Pilatus: "Quod scripsi, scripsi."
%
% Ioannes 19:21-22

