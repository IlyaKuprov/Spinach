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

% Isotope specification
sys.isotopes={};
for n=1:numel(sys_parts)
    sys.isotopes=[sys.isotopes sys_parts{n}.isotopes];
    sys_parts{n}=rmfield(sys_parts{n},'isotopes');
end

% Labels
sys.labels={};
for n=1:numel(sys_parts)
    if isfield(sys_parts{n},'labels')
        sys.labels=[sys.labels sys_parts{n}.labels];
        sys_parts{n}=rmfield(sys_parts{n},'labels');
    end 
end
if isempty(sys.labels), sys=rmfield(sys,'labels'); end

% Cartesian coordinates
inter.coordinates={};
for n=1:numel(inter_parts)
    inter.coordinates=[inter.coordinates; inter_parts{n}.coordinates];
    inter_parts{n}=rmfield(inter_parts{n},'coordinates');
end

% Zeeman interactions
inter.zeeman.matrix={};
for n=1:numel(inter_parts)
    inter.zeeman.matrix=[inter.zeeman.matrix inter_parts{n}.zeeman.matrix];
    inter_parts{n}.zeeman=rmfield(inter_parts{n}.zeeman,'matrix');
    if isempty(fieldnames(inter_parts{n}.zeeman))
        inter_parts{n}=rmfield(inter_parts{n},'zeeman');
    end
end

% Scalar couplings
inter.coupling.scalar=[];
for n=1:numel(inter_parts)
    inter.coupling.scalar=blkdiag(inter.coupling.scalar,...
                                  cell2mat(inter_parts{n}.coupling.scalar));
    inter_parts{n}.coupling=rmfield(inter_parts{n}.coupling,'scalar');
    if isempty(fieldnames(inter_parts{n}.coupling))
        inter_parts{n}=rmfield(inter_parts{n},'coupling');
    end
end
inter.coupling.scalar=num2cell(inter.coupling.scalar);

% Catch leftover subfields
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

end

% Consistency enforcement
function grumble(sys_parts,inter_parts)
if (~iscell(sys_parts))||(~iscell(inter_parts))
    error('both inputs must be cell arrays of data structures.');
end
end

% Dicebant ergo Pilato pontifices Judaeorum: "Noli scribere
% 'Rex Judaeorum' sed quia ipse dixit 'Rex sum Judaeorum'."
% Respondit Pilatus: "Quod scripsi, scripsi."
%
% Ioannes 19:21-22

