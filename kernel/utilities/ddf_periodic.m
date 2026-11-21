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
if (~isa(npts,'double'))||(~isequal(size(npts),[1 3]))||(~isreal(npts))||...
   (any(~isfinite(npts(:))))||(any(npts(:)~=fix(npts(:))))||(any(npts(:)<=0))
    error('npts must be a 1-by-3 double row of positive integers.');
end
if (~isa(spacing,'double'))||(~isequal(size(spacing),[1 3]))||(~isreal(spacing))||...
   (any(~isfinite(spacing(:))))||(any(spacing(:)<=0))
    error('spacing must be a real finite positive double 1-by-3 array.');
end
if (~isa(demag_factor,'double'))||(~isscalar(demag_factor))||(~isreal(demag_factor))||...
   (demag_factor<0)||(demag_factor>1)||(~isfinite(demag_factor))
    error('demag_factor must be a real double scalar in the interval [0,1].');
end
end


