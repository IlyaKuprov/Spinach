% Validates an explicitly supplied spherical integration grid. Syntax:
%
%                         grid_check(sph_grid)
%
% Parameters:
%
%     sph_grid - scalar structure with alphas, betas, gammas, and
%                weights fields; all four must be non-empty real
%                finite numeric vectors of equal length. Angles
%                are in radians; weights are used as supplied,
%                without normalisation or a positivity constraint.
%
% ilya.kuprov@weizmann.ac.il

function grid_check(sph_grid)

% Check the user-supplied quadrature data
grumble(sph_grid);

end

% Consistency enforcement
function grumble(sph_grid)
if (~isstruct(sph_grid))||(~isscalar(sph_grid))
    error('spherical grid must be a scalar structure.');
end
fields={'alphas','betas','gammas','weights'};
for n=1:numel(fields)
    if ~isfield(sph_grid,fields{n})
        error('spherical grid must contain alphas, betas, gammas, and weights.');
    end
    values=sph_grid.(fields{n});
    if (~isnumeric(values))||(~isreal(values))||(~isvector(values))||...
       isempty(values)||any(~isfinite(values))
        error('spherical grid fields must be non-empty real finite numeric vectors.');
    end
    if numel(values)~=numel(sph_grid.alphas)
        error('spherical grid fields must contain the same number of elements.');
    end
end
end


