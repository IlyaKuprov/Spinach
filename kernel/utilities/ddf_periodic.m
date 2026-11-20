% Secular dipolar field operator for a periodic Cartesian grid. Syntax:
%
%              symbol=ddf_periodic(npts,spacing,demag_factor)
%
% Parameters:
%
%    npts         - 1-by-3 Cartesian voxel counts in ndgrid order
%    spacing      - 1-by-3 voxel edge lengths in metres
%    demag_factor - axial demagnetising factor for the uniform mode
%
% Outputs:
%
%    symbol       - longitudinal Fourier multiplier for ddf_cartesian()

function symbol=ddf_periodic(npts,spacing,demag_factor)

% Check consistency
grumble(npts,spacing,demag_factor);

% Assemble reciprocal axes in native FFT order
axes_k=cell(1,3);
for n=1:3
    modes=-floor(npts(n)/2):(ceil(npts(n)/2)-1);
    axes_k{n}=ifftshift(modes)/(npts(n)*spacing(n));
end
[kx,ky,kz]=ndgrid(axes_k{:});

% Separate the uniform boundary condition from the secular angular tensor
norm_sq=kx.^2+ky.^2+kz.^2;
norm_sq(1)=1;
symbol=1/3-kz.^2./norm_sq;
symbol(1)=1/3-demag_factor;

end

% Consistency enforcement
function grumble(npts,spacing,demag_factor)
validateattributes(npts,{'double'},{'size',[1 3],'real','finite','integer','positive'});
validateattributes(spacing,{'double'},{'size',[1 3],'real','finite','positive'});
validateattributes(demag_factor,{'double'},{'scalar','real','>=',0,'<=',1});
end


