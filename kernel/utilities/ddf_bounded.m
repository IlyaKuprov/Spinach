% Secular dipolar field operator for isolated rectangular voxels. Syntax:
%
%                   symbol=ddf_bounded(npts,spacing)
%
% Parameters:
%
%    npts    - 1-by-3 voxel counts in Cartesian ndgrid order
%    spacing - 1-by-3 voxel edge lengths in metres
%
% Outputs:
%
%    symbol  - longitudinal Fourier multiplier for ddf_cartesian()
%
% Sources are piecewise constant and targets are cell centres. The
% opposite faces normal to Z subtend signed solid angles there; the
% self-cell includes the Lorentz spherical cavity correction (+1/3).

function symbol=ddf_bounded(npts,spacing)

% Check consistency
grumble(npts,spacing);

% Embed all source-target offsets without periodic copies
offsets=cell(1,3);
for n=1:3
    offsets{n}=[0:(npts(n)-1) (1-npts(n)):-1]*spacing(n);
end
[x,y,z]=ndgrid(offsets{:});

% Integrate the two face charges by their signed corner solid angles
weights=zeros(size(x));
for a=[-1 1]
    for b=[-1 1]
        for c=[-1 1]
            dx=x+a*spacing(1)/2;
            dy=y+b*spacing(2)/2;
            dz=z+c*spacing(3)/2;
            radius=sqrt(dx.^2+dy.^2+dz.^2);
            weights=weights-a*b*c*atan((dx.*dy)./(dz.*radius))/(4*pi);
        end
    end
end

% Remove the continuum contact field and prepare the linear convolution
weights(1)=weights(1)+1/3;
symbol=real(fftn(weights));

end

% Consistency enforcement
function grumble(npts,spacing)
if (~isa(npts,'double'))||(~isequal(size(npts),[1 3]))||(~isreal(npts))||...
   (any(~isfinite(npts(:))))||(any(npts(:)~=fix(npts(:))))||(any(npts(:)<=0))
    error('npts must be a 1-by-3 double row of positive integers.');
end
if (~isa(spacing,'double'))||(~isequal(size(spacing),[1 3]))||(~isreal(spacing))||...
   (any(~isfinite(spacing(:))))||(any(spacing(:)<=0))
    error('spacing must be a real finite positive double 1-by-3 array.');
end
end


