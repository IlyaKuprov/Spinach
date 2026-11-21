% Nine-ring Ci hemisphere quadrature from doi:10.5194/mr-6-15-2025.
% The third row of the active ZYZ rotation is the molecular-frame field.
%
% Syntax:
%
%                        grid=trityl_grid()
%
% Outputs:
%
%   grid - column vectors of active ZYZ Euler angles (rad) and weights
%
% ilya.kuprov@weizmann.ac.il

function grid=trityl_grid()

% Integrate the cells of the nine-ring hemisphere quadrature
knots=9;
step=pi/(2*(knots-1));
grid.alphas=0;
grid.betas=0;
grid.gammas=0;
grid.weights=1-cos(step/2);
for n=1:knots-1

    % Use four additional azimuths for each step away from the pole
    theta=n*step;
    azimuth=2*pi*(0:4*n-1).'/(4*n);
    grid.alphas=[grid.alphas; zeros(4*n,1)];
    grid.betas=[grid.betas; theta*ones(4*n,1)];
    grid.gammas=[grid.gammas; pi-azimuth];
    lower=theta-step/2;
    upper=min(theta+step/2,pi/2);
    grid.weights=[grid.weights; repmat((cos(lower)-cos(upper))/(4*n),4*n,1)];

end

end


