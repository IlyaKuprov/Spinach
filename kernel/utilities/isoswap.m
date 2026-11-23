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
% Note: nuclear quadrupolar couplings are rescaled using the quadrupole
%       moments and the spin quantum numbers of the two isotopes. Other
%       quadratic couplings and higher order couplings are wiped and a
%       warning is printed - those are not transferable.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=isoswap.m>

function [sys,inter]=isoswap(sys,inter,spins,new_iso)

% Grumbler missing
grumble(sys,inter,spins,new_iso);

% Wipe high-rank couplings
if isfield(inter,'giant')
    for n=spins(:)'
        if ~isempty(inter.giant.coeff{n})
            inter.giant.coeff{n}={}; inter.giant.euler{n}={};
            disp(['WARNING: high-rank coupling for spin ' int2str(n) ...
                  ' wiped (not transferable).']);
        end
    end
end

% Loop over spins
for n=spins(:)'

    % Get the gamma ratio
    gamma_ratio=spin(new_iso)/spin(sys.isotopes{n});

    % Get the ratio of quadrupolar coupling prefactors, Q/(2I(2I-1))
    [~,old_mult,old_data]=spin(sys.isotopes{n}); [~,new_mult,new_data]=spin(new_iso);
    quad_ratio=(new_data.quad_moment/((new_mult-1)*(new_mult-2)))/...
               (old_data.quad_moment/((old_mult-1)*(old_mult-2)));
    transferable=isscalar(quad_ratio)&&isfinite(quad_ratio)&&(quad_ratio~=0);

    % Scale or wipe quadratic couplings - eigensystem representation
    if isfield(inter,'coupling')&&isfield(inter.coupling,'eigs')&&(~isempty(inter.coupling.eigs{n,n}))
        if transferable
            inter.coupling.eigs{n,n}=quad_ratio*inter.coupling.eigs{n,n};
        else
            inter.coupling.eigs{n,n}=[]; inter.coupling.euler{n,n}=[];
            disp(['WARNING: quadratic coupling for spin ' int2str(n) ...
                  ' wiped (not transferable).']);
        end
    end

    % Scale or wipe quadratic couplings - matrix representation
    if isfield(inter,'coupling')&&isfield(inter.coupling,'matrix')&&(~isempty(inter.coupling.matrix{n,n}))
        if transferable
            inter.coupling.matrix{n,n}=quad_ratio*inter.coupling.matrix{n,n};
        else
            inter.coupling.matrix{n,n}=[];
            disp(['WARNING: quadratic coupling for spin ' int2str(n) ...
                  ' wiped (not transferable).']);
        end
    end

    % Loop over coupling partners
    for k=setdiff(1:numel(sys.isotopes),n)

        % Scale eigensystem representation
        if isfield(inter,'coupling')&&isfield(inter.coupling,'eigs')
            inter.coupling.eigs{n,k}=gamma_ratio*inter.coupling.eigs{n,k};
            inter.coupling.eigs{k,n}=gamma_ratio*inter.coupling.eigs{k,n};
        end

        % Scale matrix representation
        if isfield(inter,'coupling')&&isfield(inter.coupling,'matrix')
            inter.coupling.matrix{n,k}=gamma_ratio*inter.coupling.matrix{n,k};
            inter.coupling.matrix{k,n}=gamma_ratio*inter.coupling.matrix{k,n};
        end

        % Scale scalar couplings
        if isfield(inter,'coupling')&&isfield(inter.coupling,'scalar')
            inter.coupling.scalar{n,k}=gamma_ratio*inter.coupling.scalar{n,k};
            inter.coupling.scalar{k,n}=gamma_ratio*inter.coupling.scalar{k,n};
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
    if any(cellfun(@spin,sys.isotopes(spins))==0)
        error('interactions of zero-gamma particles cannot be rescaled.');
    end
end

% Он из полка был изгнан за дуэль
% Или за то, что не был на дуэли.
% 
% Михаил Лермонтов

