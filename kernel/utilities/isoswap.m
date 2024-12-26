% Makes isotope replacements in the input structures. All interactions
% are automatically scaled as necessary. Syntax:
%
%             [sys,inter]=isoswap(sys,inter,spins,new_iso)
%
% Parameters:
%
%   sys, inter  - Spinach input data structures
%
%   spins       - a vector of integers specifying
%                 spin numbers
%
%   new_iso     - character string specifying 
%                 the new isotope
%
% Output:
%
%   sys, inter  - Spinach input data structures
%
% Note: quadratic and higher order couplings are wiped and a warning
%       is printed - those are not transferable.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=isoswap.m>

function [sys,inter]=isoswap(sys,inter,spins,new_iso)

% Grumbler missing
grumble(sys,inter,spins,new_iso);

% Wipe quadratic couplings - eigensystem representation
if isfield(inter,'coupling')&&isfield(inter.coupling,'eigs')
    for n=spins(:)'
        if ~isempty(inter.coupling.eigs{n,n})
            inter.coupling.eigs{n,n}=[]; inter.coupling.euler{n,n}=[];
            disp(['WARNING: quadratic coupling for spin ' int2str(n) ...
                  ' wiped (not transferable).']);
        end
    end
end

% Wipe quadratic couplings - matrix representation
if isfield(inter,'coupling')&&isfield(inter.coupling,'matrix')
    for n=spins(:)'
        if ~isempty(inter.coupling.matrix{n,n})
            inter.coupling.matrix{n,n}=[];
            disp(['WARNING: quadratic coupling for spin ' int2str(n) ...
                  ' wiped (not transferable).']);
        end
    end
end

% Wipe high-rank couplings
if isfield(inter,'giant')
    for n=spins(:)'
        if ~isempty(inter.giant.coef{n})
            inter.giant.coef{n}={}; inter.giant.euler{n}={};
            disp(['WARNING: high-rank coupling for spin ' int2str(n) ...
                  ' wiped (not transferable).']);
        end
    end
end

% Loop over spins
for n=spins(:)'

    % Get the gamma ratio
    gamma_ratio=spin(new_iso)/spin(sys.isotopes{n});

    % Loop over coupling partners
    for k=setdiff(1:numel(sys.isotopes),n)

        % Scale eigensystem representation
        if isfield(inter.coupling,'eigs')
            inter.coupling.eigs{n,k}=gamma_ratio*inter.coupling.eigs{n,k};
            inter.coupling.eigs{k,n}=gamma_ratio*inter.coupling.eigs{k,n};
        end

        % Scale matrix representation
        if isfield(inter.coupling,'matrix')
            inter.coupling.matrix{n,k}=gamma_ratio*inter.coupling.matrix{n,k};
            inter.coupling.matrix{k,n}=gamma_ratio*inter.coupling.matrix{k,n};
        end

    end

    % Replace the isotope string
    sys.isotopes{n}=new_iso;

end

end

% Consistency enforcement
function grumble(sys,inter,spins,new_iso)
    if ~isfield(sys,'isotopes')
        error('isotope information missing from sys structure.');
    end
    if ~isstruct(inter)
        error('inter must be a Spinach interactions data structure.');
    end
    if (~isnumeric(spins))||(~isreal(spins))||(~isvector(spins))||...
       (any(mod(spins,1)~=0,'all'))||(any(spins>numel(sys.isotopes),'all'))
        error('spins must be a vector of integer spin indices.');
    end
    if ~ischar(new_iso)
        error('new_iso must be a character string.');
    end
end

% Он из полка был изгнан за дуэль
% Или за то, что не был на дуэли.
% 
% Михаил Лермонтов

