% Generates repulsion grids on a unit hypersphere. See the paper by
% Bak and Nielsen (http://dx.doi.org/10.1006/jmre.1996.1087) to get
% further information on the algorithm involved. Syntax:
%
%    [alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)
%
% Parameters:
%
%     npoints - number of points in the resulting spherical grid
%
%       ndims - hypersphere dimension: 2 returns a single-angle 
%               (beta) grid, 3 returns a two-angle grid (alpha,
%               beta), 4 returns a three-angle (alpha,beta,gam-
%               ma) spherical grid
%
%       niter - number of repulsion interations (simple clipped
%               gradient descent at the moment)
%
% Outputs:
%
%      alphas - alpha Euler angles of the grid, in radians,
%               zeros for two-angle grids
%
%       betas - beta Euler angles of the grid, in radians
%
%      gammas - gamma Euler angles of the grid, in radians,
%               zeros for single-angle grids
%
%     weights - point weights of the grid
% 
% Note: uniform weights are assigned at the moment, use the supp-
%       lied SHREWD function to generate optimal weights.
%
% i.kuprov@soton.ac.uk
% fmentink@magnet.fsu.edu
%
% <https://spindynamics.org/wiki/index.php?title=repulsion.m>

function [alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)

% Check consistency
grumble(npoints,ndims,niter);

% Generate guess points
R=rand(ndims,npoints)-0.5;

% Start the repulsion loop
for m=1:niter
    
    % Make space for distance vector array
    Rd=reshape(R,[ndims 1 npoints]);
    
    % Compute distance vector array
    dist_vecs=permute(Rd,[1 2 3])-...
              permute(Rd,[1 3 2]);
    
    % Normalise distance vectors
    dist_vecs=dist_vecs./sqrt(sum(dist_vecs.^2,1));
    
    % Eliminate self-interactions
    dist_vecs(~isfinite(dist_vecs))=0;
    
    % Get scalar products
    scalar_prods=reshape(R'*R,[1 npoints npoints]);
    
    % Get tangent forces
    F=sum(dist_vecs.*scalar_prods,3);
    
    % Move under tangent forces
    R_new=R-ndims*F/npoints;
    
    % Reproject onto unit sphere
    R_new=R_new./sqrt(sum(R_new.^2,1));
    
    % Report the difference
    max_diff=max(sqrt(sum((R-R_new).^2,2)));
    disp(['Iteration ' num2str(m) ...
          ', maximum point displacement: ' num2str(max_diff)]);
    
    % Close the loop
    R=R_new;
    
end

% Get points
switch ndims
    
    case 2
        
        % In 2D case return polar angles
        [phi,~]=cart2pol(R(1,:),R(2,:));
        betas=phi'; alphas=0*betas; gammas=0*betas;
        
    case 3
        
        % In 3D case return spherical angles
        [phi,theta,~]=cart2sph(R(1,:),R(2,:),R(3,:));
        betas=theta'+pi/2; gammas=phi'; alphas=0*gammas;
        
        % Display the grid
        if nargout==0
            figure(); plot3(cos(gammas).*sin(betas),...
                            sin(gammas).*sin(betas),...
                            cos(betas),'r.');
            axis square; axis tight; kgrid; box on;
        end
    
    case 4
        
        % In 4D case return Euler angles
        [alphas,betas,gammas]=quat2angle(R','ZYZ');
       
end

% Get weights
weights=ones(npoints,1)/npoints;

end

% Consistency enforcement
function grumble(npoints,ndims,niter)
if (~isnumeric(npoints))||(~isreal(npoints))||(numel(npoints)~=1)||(npoints<1)||(mod(npoints,1)~=0)
    error('npoints must be a positive integer.');
end
if (~isnumeric(niter))||(~isreal(niter))||(numel(niter)~=1)||(niter<1)||(mod(niter,1)~=0)
    error('niter must be a positive integer.');
end
if (~isnumeric(ndims))||(~isreal(ndims))||(numel(ndims)~=1)||(~ismember(ndims,[2 3 4]))
    error('ndims must be 2, 3 or 4.');
end
end
                                                          
% Let us beware of saying that there are laws in nature. There 
% are only necessities: there is nobody who commands, nobody
% who obeys, nobody who trespasses.
%
% Friedrich Nietzsche

