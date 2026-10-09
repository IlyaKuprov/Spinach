% Rigid-molecule rotational diffusion estimate from an ELM-inspired
% hydrated surface ellipsoid. Syntax:
%
%   [tau,D,axes_len]=rotcorr(xyz,radii,layer,probe,temp,visc,npoints)
%
% Parameters:
%
%   xyz       - real finite floating-point N-by-3 Cartesian
%               coordinates, Angstrom
%   radii     - positive finite floating-point N-by-1 atomic radii,
%               Angstrom;
%               choose radii appropriate to the atom list and model
%   layer     - nonnegative hydration-shell thickness, Angstrom
%   probe     - nonnegative solvent probe radius, Angstrom
%   temp      - positive absolute temperature, Kelvin
%   visc      - positive dynamic solvent viscosity, Pa s
%   npoints   - integer number of directions per atom, at least 6;
%               increase until the requested outputs converge
%
% Outputs:
%
%   tau       - isotropic-equivalent rank-2 time, 1/(2*trace(D)), s
%   D         - 3-by-3 rotational diffusion tensor in the input
%               Cartesian frame, inverse seconds
%   axes_len  - 3-by-1 ellipsoid semi-axes, descending, Angstrom
%
% The model retains exposed contact points on spheres of radius
% radii+layer, tested with a rolling solvent probe. Equal-solid-angle
% samples carry sphere-area weights. Surface covariance gives the
% semi-axes sqrt(3*eigenvalues); Perrin integrals give diffusion.
% This is a sampled-contact-surface approximation, not a reproduction
% of SURF triangulation, and excludes re-entrant surface patches.
% Numerical cost is O(N^2+N*npoints*neighbours), with O(N+npoints)
% workspace. The inputs must describe one rigid body in dilute fluid.
%
% Based on Ryabov et al., JACS 128 (2006), 15432-15444,
% https://doi.org/10.1021/ja062715t. Hydration and radii are model
% parameters, not quantities determined by coordinates. Ellipsoidal
% reduction loses detailed shape; flexible or strongly non-globular
% molecules require a more complete hydrodynamic treatment. tau is
% not a universal single-exponential correlation time for anisotropic
% rotation. D retains the anisotropy needed for that treatment.
%
% talos@spindynamics.org
%
% <https://spindynamics.org/wiki/index.php?title=rotcorr.m>

function [tau,D,axes_len]=rotcorr(xyz,radii,layer,probe,temp,visc,npoints)

% Check consistency
grumble(xyz,radii,layer,probe,temp,visc,npoints);

% Remove an irrelevant translation before accumulating moments
xyz=double(xyz); radii=double(radii);
layer=double(layer); probe=double(probe);
temp=double(temp); visc=double(visc); npoints=double(npoints);
xyz=xyz-mean(xyz,1); shell=radii+layer;

% Build an equal-solid-angle Fibonacci sphere
z=1-2*((1:npoints)'-0.5)/npoints;
phi=pi*(3-sqrt(5))*(0:npoints-1)';
dirs=[sqrt(1-z.^2).*cos(phi) sqrt(1-z.^2).*sin(phi) z];

% Accumulate area-weighted exposed surface moments
area=0; first=zeros(1,3); second=zeros(3);
for n=1:size(xyz,1)

    % Only overlapping expanded spheres can occlude this atom
    delta=xyz-xyz(n,:); dist=sqrt(sum(delta.^2,2));
    neighbours=find((dist<shell+shell(n)+2*probe)&...
                    ((1:size(xyz,1))'~=n));
    centres=xyz(n,:)+(shell(n)+probe)*dirs;
    exposed=true(npoints,1);
    for k=neighbours'
        delta=centres-xyz(k,:);
        exposed=exposed&(sum(delta.^2,2)>=(shell(k)+probe)^2);
    end

    % Integrate contact-surface moments with sphere-area weights
    points=xyz(n,:)+shell(n)*dirs(exposed,:);
    weight=shell(n)^2/npoints;
    area=area+weight*size(points,1);
    first=first+weight*sum(points,1);
    second=second+weight*(points'*points);
end

% Obtain the surface-covariance ellipsoid
if area==0
    error('no exposed surface samples; increase npoints.');
end
centre=first/area; covar=second/area-centre'*centre;
[frame,values]=eig(covar,'vector');
[values,order]=sort(values,'descend'); frame=frame(:,order);
if any(values<=0)
    error('surface covariance is singular; increase npoints.');
end
axes_len=sqrt(3*values);

% Evaluate dimensionless Perrin integrals to avoid SI underflow
scale=max(axes_len); squares=(axes_len/scale).^2;
ints=zeros(3,1);
for n=1:3
    ints(n)=integral(@(u)1./((squares(n)+u).*...
        sqrt((squares(1)+u).*(squares(2)+u).*(squares(3)+u))),...
        0,Inf,'RelTol',1e-10,'AbsTol',0);
end

% Convert ellipsoid friction to diffusion in the original frame
terms=squares.*ints;
rates=3*1.380649e-23*temp/(16*pi*visc*(scale*1e-10)^3)*...
      (sum(terms)-terms)./(sum(squares)-squares);
D=frame*diag(rates)*frame'; tau=1/(2*sum(rates));

end

% Validate the molecular, solvent, and quadrature specifications
function grumble(xyz,radii,layer,probe,temp,visc,npoints)
if ~isfloat(xyz)||~isreal(xyz)||size(xyz,2)~=3||...
   ~ismatrix(xyz)||isempty(xyz)||any(~isfinite(xyz),'all')
    error('xyz must be a nonempty real finite N-by-3 array.');
end
if ~isfloat(radii)||~isreal(radii)||...
   ~isequal(size(radii),[size(xyz,1) 1])||...
   any(~isfinite(radii))||any(radii<=0)
    error('radii must be a positive finite N-by-1 array.');
end
if ~isnumeric(layer)||~isreal(layer)||~isscalar(layer)||...
   ~isfinite(layer)||layer<0
    error('layer must be a nonnegative finite scalar.');
end
if ~isnumeric(probe)||~isreal(probe)||~isscalar(probe)||...
   ~isfinite(probe)||probe<0
    error('probe must be a nonnegative finite scalar.');
end
if ~isnumeric(temp)||~isreal(temp)||~isscalar(temp)||...
   ~isfinite(temp)||temp<=0
    error('temp must be a positive finite scalar.');
end
if ~isnumeric(visc)||~isreal(visc)||~isscalar(visc)||...
   ~isfinite(visc)||visc<=0
    error('visc must be a positive finite scalar.');
end
if ~isnumeric(npoints)||~isreal(npoints)||~isscalar(npoints)||...
   ~isfinite(npoints)||npoints<6||mod(npoints,1)~=0
    error('npoints must be an integer of at least 6.');
end
end


