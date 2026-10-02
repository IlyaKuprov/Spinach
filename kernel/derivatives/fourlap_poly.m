% Switch-aware periodic Fourier Laplacian. Syntax:
%
%                  L=fourlap_poly(spin_system,npoints,extents)
%
% Parameters:
%
%    spin_system - Spinach system with the sys.enable option list
%    npoints     - row vector of grid sizes in one to three dimensions
%    extents     - row vector of positive periods in matching order
%
% Outputs:
%
%    L           - Fourier Laplacian on vectorised X, Y, Z grid data
%
% With polyadics enabled, each second derivative uses inverse FFT,
% a numeric CPU multiplier, and FFT factors. Otherwise, fourlap
% supplies its existing explicit matrix. The two-argument fourlap
% API is unchanged. Periodic boundary conditions are required.
% Upload an implicit result once with gpuArray for GPU actions.
%
% ilya.kuprov@weizmann.ac.il

function L=fourlap_poly(spin_system,npoints,extents)

% Check consistency
grumble(spin_system,npoints,extents);

% Retain the existing explicit matrix when polyadics are disabled
if ~ismember('polyadic',spin_system.sys.enable)
    L=fourlap(npoints,extents); return;
end

% Build a sum of second derivatives in reverse tensor order
terms=cell(1,numel(npoints));
for axis=1:numel(npoints)
    cores=cell(1,numel(npoints));
    for dim=1:numel(npoints)
        if dim==axis
            [~,cores{dim}]=fourdif(spin_system,npoints(dim),2);
            cores{dim}=(2*pi/extents(dim))^2*cores{dim};
        else
            cores{dim}=opium(npoints(dim),1);
        end
    end
    terms{axis}=fliplr(cores);
end
L=polyadic(terms);

end

% Consistency enforcement
function grumble(spin_system,npoints,extents)
if (~isstruct(spin_system))||(~isfield(spin_system,'sys'))||...
   (~isfield(spin_system.sys,'enable'))||(~iscell(spin_system.sys.enable))
    error('spin_system.sys.enable must be a cell array of option names.');
end
if (~isnumeric(npoints))||(~isreal(npoints))||(~isrow(npoints))||...
   (numel(npoints)<1)||(numel(npoints)>3)||any(~isfinite(npoints))||...
   any(npoints<1)||any(mod(npoints,1)~=0)
    error('npoints must be a row of one to three finite positive integers.');
end
if (~isnumeric(extents))||(~isreal(extents))||(~isrow(extents))||...
   (numel(extents)~=numel(npoints))||any(~isfinite(extents))||any(extents<=0)
    error('extents must be a matching row of finite positive real periods.');
end
end


