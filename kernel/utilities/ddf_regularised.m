% Isolated softened dipolar kernel for midpoint quadrature. Syntax:
%
%           symbol=ddf_regularised(npts,spacing,radius,self_term)
%
% Parameters:
%
%    npts      - 1-by-3 Cartesian voxel counts
%    spacing   - 1-by-3 voxel side lengths, metres
%    radius    - positive softening length, metres
%    self_term - 'omit' (zero diagonal) or 'regularised' (retain it)
%
% Outputs:
%
%    symbol    - dimensionless longitudinal multiplier for ddf_cartesian
%
% At coincidence the squared direction cosine is defined as zero.
% This comparison kernel is not the cell-integrated Lorentz cavity.
% To obtain the BlochDDFGalerkin angular-frequency field for scalar
% k_ddf, multiply the output of ddf_cartesian by -4*pi*k_ddf/mu0.

function symbol=ddf_regularised(npts,spacing,radius,self_term)

% Check consistency
grumble(npts,spacing,radius,self_term);

% Construct all distinct displacements for a linear convolution
offsets=cell(1,3);
for n=1:3
    offsets{n}=[0:(npts(n)-1) (1-npts(n)):-1]*spacing(n);
end
[x,y,z]=ndgrid(offsets{:});
squared=x.^2+y.^2+z.^2;

% Form the softened axial tensor with midpoint source volumes
cos_sq=zeros(size(z));
away=(squared>0);
cos_sq(away)=z(away).^2./squared(away);
weights=prod(spacing)*(3*cos_sq-1)./(4*pi*(squared+radius^2).^(3/2));
if strcmp(self_term,'omit')
    weights(1)=0;
end
symbol=real(fftn(weights));

end

% Consistency enforcement
function grumble(npts,spacing,radius,self_term)
validateattributes(npts,{'double'},{'size',[1 3],'real','finite','integer','positive'});
validateattributes(spacing,{'double'},{'size',[1 3],'real','finite','positive'});
validateattributes(radius,{'double'},{'scalar','real','finite','positive'});
if ~ischar(self_term)||~ismember(self_term,{'omit','regularised'})
    error('self_term must be omit or regularised.');
end
end


