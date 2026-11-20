% Face-flux diffusion rate operator on a reflecting voxel mask. Syntax:
%
%               rate=ddf_diffusion(mask,spacing,diffusivity)
%
% Parameters:
%
%    mask        - logical Cartesian array marking permeable voxels
%    spacing     - 1-by-3 voxel edge lengths in metres
%    diffusivity - uniform isotropic diffusivity in m^2/s
%
% Outputs:
%
%    rate        - sparse numel(mask)-square rate matrix in 1/s,
%                  with zero flux across closed faces and grid edges

function rate=ddf_diffusion(mask,spacing,diffusivity)

% Check consistency
grumble(mask,spacing,diffusivity);

% Assemble undirected face conductances without crossing a mask boundary
count=numel(mask);
indices=reshape(1:count,size(mask));
rate=sparse(count,count);
for n=1:3
    left=repmat({':'},1,3);
    right=left;
    left{n}=1:(size(mask,n)-1);
    right{n}=2:size(mask,n);
    active=mask(left{:})&mask(right{:});
    a=indices(left{:});
    b=indices(right{:});
    a=a(active);
    b=b(active);
    faces=sparse(a(:),b(:),diffusivity/spacing(n)^2,count,count);
    rate=rate+faces+faces.';
end

% Balance each voxel's incoming and outgoing fluxes
rate=rate-spdiags(sum(rate,2),0,count,count);

end

% Consistency enforcement
function grumble(mask,spacing,diffusivity)
if ~islogical(mask)||isempty(mask)||(ndims(mask)>3)
    error('mask must be a nonempty logical Cartesian array.');
end
if (~isa(spacing,'double'))||(~isequal(size(spacing),[1 3]))||(~isreal(spacing))||...
   (any(~isfinite(spacing(:))))||(any(spacing(:)<=0))
    error('spacing must be a real finite positive double 1-by-3 array.');
end
if (~isa(diffusivity,'double'))||(~isscalar(diffusivity))||(~isreal(diffusivity))||...
   (~isfinite(diffusivity))||(diffusivity<0)
    error('diffusivity must be a real finite nonnegative double scalar.');
end
end


