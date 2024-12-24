% Harmonic oscillator infrastructure in 1D. Syntax:
%
%           [H_oscl,X_oscl,xgrid]=oscillator(parameters)
%
% Parameters:
%
%   parameters.frc_cnst    - force constant, N/m
%
%   parameters.par_mass    - particle mass, kg
%
%   parameters.grv_cnst    - gravitational acceleration. m/s^2
%
%   parameters.n_points    - number of discretization points
%
%   parameters.box_size    - oscillator box size, m
%
% Outputs:
%
%   H_oscl   -   oscillator Hamiltonian, Joules
%
%   X_oscl   -   oscillator X operator, m
%
%   xgrid    -   X coordinate grid, m
%
% Note: gravitation is directed along the X axis. Finite difference
%       derivative operators are used.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=oscillator.m>

function [H_oscl,X_oscl,xgrid]=oscillator(parameters)

% Check consistency
grumble(parameters);

% Second derivative operator
d2_dx2=((parameters.n_points-1)/parameters.box_size)^2*fdmat(parameters.n_points,5,2);

% Coordinate grid
xgrid=linspace(-parameters.box_size/2,parameters.box_size/2,parameters.n_points)';

% Coordinate operator
X_oscl=spdiags(xgrid,0,parameters.n_points,parameters.n_points);

% Boundary conditions
d2_dx2(:,1)=0; d2_dx2(:,end)=0; d2_dx2(1,:)=0; d2_dx2(end,:)=0;

% Hamiltonian
H_oscl=-(1/(2*parameters.par_mass))*d2_dx2+...
        (parameters.frc_cnst/2)*X_oscl^2+...
         parameters.par_mass*parameters.grv_cnst*X_oscl;

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'frc_cnst')
    error('force constant must be specified in parameters.frc_cnst variable.');
end
if (~isnumeric(parameters.frc_cnst))||(~isreal(parameters.frc_cnst))||...
   (numel(parameters.frc_cnst)~=1)||(parameters.frc_cnst<=0)
    error('parameters.frc_cnst must be a positive real number.');
end
if ~isfield(parameters,'par_mass')
    error('particle mass must be specified in parameters.par_mass variable.');
end
if (~isnumeric(parameters.par_mass))||(~isreal(parameters.par_mass))||...
   (numel(parameters.par_mass)~=1)||(parameters.par_mass<=0)
    error('parameters.par_mass must be a positive real number.');
end
if ~isfield(parameters,'grv_cnst')
    error('gravitational constant must be specified in parameters.grv_cnst variable.');
end
if (~isnumeric(parameters.grv_cnst))||(~isreal(parameters.grv_cnst))||...
   (numel(parameters.grv_cnst)~=1)
    error('parameters.grv_cnst must be a real number.');
end
if ~isfield(parameters,'n_points')
    error('number of grid points must be specified in parameters.n_points variable.');
end
if (~isnumeric(parameters.n_points))||(~isreal(parameters.n_points))||...
   (numel(parameters.n_points)~=1)||mod(parameters.n_points,1)||(parameters.n_points<1)
    error('parameters.n_points must be a positive real integer.');
end
if ~isfield(parameters,'box_size')
    error('box size must be specified in parameters.box_size variable.');
end
if (~isnumeric(parameters.box_size))||(~isreal(parameters.box_size))||...
   (numel(parameters.box_size)~=1)||(parameters.box_size<=0)
    error('parameters.box_size must be a positive real number.');
end
end

% Just in terms of allocation of time resources, religion is
% not very efficient. There's a lot more I could be doing on
% a Sunday morning.
%
% Bill Gates

