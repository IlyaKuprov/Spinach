% Plot orientation-selection weights on a spherical grid. Syntax:
%
%             plot_orisel(sph_grid,selection_weights,transition_idx)
%
% Parameters:
%
%   sph_grid          - alphas, betas, gammas, weights: grid columns
%   selection_weights - nrows-by-norient non-negative intensity weights
%   transition_idx    - one row of selection_weights to plot
%
% The chosen row is scaled to unit maximum; an all-zero row is shown
% at the bottom of the colour scale. No quadrature weights are applied.
% The sphere direction is (sin(beta)cos(gamma),sin(beta)sin(gamma),
% cos(beta)). A full sphere with distinct directions is required:
% grids that repeat beta/gamma for multiple alpha angles are not
% supported. Existing axes are used and their hold state is preserved.
%
% alexey.bogdanov@weizmann.ac.il

function plot_orisel(sph_grid,selection_weights,transition_idx)

% Check consistency
grumble(sph_grid,selection_weights,transition_idx);

% Construct the spherical Voronoi tessellation
x=sin(sph_grid.betas).*cos(sph_grid.gammas);
y=sin(sph_grid.betas).*sin(sph_grid.gammas);
z=cos(sph_grid.betas);
[~,~,vorn]=voronoisphere([x'; y'; z']);

% Scale the requested intensity row
weights=selection_weights(transition_idx,:);
max_w=max(weights);
if max_w>0, weights=weights/max_w; end

% Draw true-colour patches and preserve the caller's hold state
cmap=parula(256); held=ishold; hold on;
for k=1:numel(vorn)
    idx=1+round(weights(k)*(size(cmap,1)-1));
    patch(vorn{k}(1,:),vorn{k}(2,:),vorn{k}(3,:),cmap(idx,:),...
          'FaceAlpha',1,'EdgeColor','none');
end
if ~held, hold off; end

% Show spatial axes and the relative-intensity colour scale
colormap(cmap); clim([0 1]); colorbar;
xlim([-1.1 1.1]); ylim([-1.1 1.1]); zlim([-1.1 1.1]);
axis equal; box on; view(3);
xticks([-1 0 1]); yticks([-1 0 1]); zticks([-1 0 1]);
xlabel('x'); ylabel('y'); zlabel('z');

end

% Check input consistency
function grumble(sph_grid,selection_weights,transition_idx)
if ~isstruct(sph_grid)
    error('sph_grid must be a spherical grid structure.');
end
fields={'alphas','betas','gammas','weights'};
for k=1:numel(fields)
    x=sph_grid.(fields{k});
    if ~isnumeric(x)||~isreal(x)||~iscolumn(x)||...
       any(~isfinite(x))||(numel(x)~=numel(sph_grid.weights))
        error('grid fields must be equal-length finite real columns.');
    end
end
if ~isnumeric(selection_weights)||~isreal(selection_weights)||...
   ~ismatrix(selection_weights)||any(~isfinite(selection_weights),'all')||...
   any(selection_weights<0,'all')||...
   (size(selection_weights,2)~=numel(sph_grid.weights))
    error('selection_weights must be finite non-negative rows matching the grid.');
end
if ~isnumeric(transition_idx)||~isreal(transition_idx)||...
   ~isscalar(transition_idx)||~isfinite(transition_idx)||...
   (mod(transition_idx,1)~=0)||(transition_idx<1)||...
   (transition_idx>size(selection_weights,1))
    error('transition_idx must identify one selection_weights row.');
end
end


