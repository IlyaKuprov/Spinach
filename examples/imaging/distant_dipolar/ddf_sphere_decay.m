% Radial diffusion in a reflecting sphere: mode projection and local
% decay for the second positive root of tan(q)=q. Matches the FD
% benchmark of J. Chem. Phys. 164, 144121 (doi:10.1063/5.0325917).
%
%            [trace,mode_error]=ddf_sphere_decay(count,radius)
%
% Inputs:
%    count      - even number of intervals across the unit-metre box
%    radius     - sphere radius in metres, strictly below 0.5
%
% Outputs:
%    trace      - columns: time (s), voxel and continuum amplitudes
%    mode_error - relative spatial L2 error at the final time

function [trace,mode_error]=ddf_sphere_decay(count,radius)

% Check consistency
grumble(count,radius);

% Keep only nodes of the unit box that can intersect the sphere
axis_vals=(-ceil(count*radius):ceil(count*radius))/count;
[x,y,z]=ndgrid(axis_vals);
radii=sqrt(x.^2+y.^2+z.^2);
mask=(radii<=radius+eps(radius));

% Initialise the radial Neumann eigenmode on the reflecting mask
root=fzero(@(q) sin(q)-q*cos(q),[2*pi 2.5*pi]);
mode=sinc((root/(pi*radius))*radii(mask));
density=0.1*mode;
rate=ddf_diffusion(mask,[1 1 1]/count,1e-3);
rate=rate(mask,mask)-speye(nnz(mask))/5;

% Match the archive step unless refinement requires explicit substeps
step=5e-4;
substeps=ceil(step*(6e-3*count^2+1/5));
times=(0:200)'*step;
trace=[times zeros(size(times)) 0.1*exp(-(1e-3*(root/radius)^2+1/5)*times)];
trace(1,2)=0.1;
for n=2:numel(times)
    for k=1:substeps
        density=density+(step/substeps)*(rate*density);
    end
    trace(n,2)=(mode'*density)/(mode'*mode);
end

% Compare the complete local magnetisation, not its radial projection alone
mode_error=norm(density-trace(end,3)*mode)/norm(trace(end,3)*mode);

end

% Consistency enforcement
function grumble(count,radius)
if (~isa(count,'double'))||(~isscalar(count))||(~isreal(count))||(~isfinite(count))||(count<=0)||...
   (count~=fix(count))
    error('count must be a real finite positive integer-valued double scalar.');
end
if (~isa(radius,'double'))||(~isscalar(radius))||(~isreal(radius))||(~isfinite(radius))||...
   (radius<=0)||(radius>=0.5)
    error('radius must be a real double scalar strictly between zero and 0.5.');
end
if mod(count,2)~=0
    error('count must be even to sample the sphere centre.');
end
end


