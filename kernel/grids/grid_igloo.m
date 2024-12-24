% Igloo grid, as per Appendix A2 of 
%
%      http://dx.doi.org/10.1016/j.jmr.2014.05.009
%
% Syntax:
%
%     [alps,bets,gams,whts,vorn]=grid_igloo(n_long)
%
% Parameters:
%
%    n_long - number of longitudes in the grid
%
% Outputs:
%
%      alps - alpha Euler angles of the grid (radians),
%             zeros because these are two-angle grids
%
%      bets - beta Euler angles of the grid (radians)
%
%      gams - gamma Euler angles of the grid (radians)
% 
%      whts - Voronoi tessellation body angle weights
%
%      vorn - a cell array of matrices containing the
%             coordinates of the vertices of the Voro-
%             noi polyhedra
%
% If no outputs are requested, a schematic is drawn.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grid_igloo.m>

function [alps,bets,gams,whts,vorn]=grid_igloo(n_long)

% Check consistency
grumble(n_long);

% Betas equispaced by longitude
bets=linspace(0,pi,n_long)';

% Preallocate three-angle blocks
alps_bets_gams=cell(n_long,3);

% Loop over longitudes
for k=1:n_long
    
    % Compute the number of lattitudes
    n_latt=floor(2*(n_long-1)*sin(bets(k))+0.5);
    n_latt=max([1 n_latt]);
    
    % Get the gamma angles of lattitudes
    gams=linspace(0,2*pi,n_latt+1)'; 
    gams(end)=[];
    
    % Fill the angle block
    alps_bets_gams{k,1}=zeros(size(gams));
    alps_bets_gams{k,2}=bets(k)*ones(size(gams));
    alps_bets_gams{k,3}=gams;
    
end

% Concatenate angle blocks
alps_bets_gams=cell2mat(alps_bets_gams);

% Extract individual angles
alps=alps_bets_gams(:,1);
bets=alps_bets_gams(:,2);
gams=alps_bets_gams(:,3);

% Weights are expensive
if (nargout>3)||(nargout==0)
    
    % Get point coordinates
    x=sin(bets).*cos(gams);
    y=sin(bets).*sin(gams);
    z=cos(bets);
    
    % Run Voronoi tessellation
    [~,~,vorn,whts]=voronoisphere([x'; y'; z']);
    
    % Get weights
    whts=whts/(4*pi);
    
end

% If no output requested, plot the grid
if nargout==0, grid_plot(x,y,z,vorn); end

end

% Consistency enforcement
function grumble(n_long)
if (~isnumeric(n_long))||(~isreal(n_long))||...
   (~isscalar(n_long))||(n_long<1)||mod(n_long,1)
    error('n_long must be a positive real integer');
end
end

% Разговор в автобусе:
% - Он не хочет, говорит что у него иссяк запал.
% - Что-что у него запало?
% - Иссяк.

