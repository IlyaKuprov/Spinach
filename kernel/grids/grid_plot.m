% Spherical quadrature grid plotter. Takes a cloud of points
% on a sphere and plots its Voronoi tessellation. Syntax:
%
%                     grid_plot(x,y,z,vorn)
%
% Parameters:
%
%   x,y,z  - column vectors containing Cartesian 
%            coordinates of grid points
%
%   vorn   - Voronoi tessellation; if this is not 
%            provided, it will be computed
%
% Outputs:
%
%   this function plots a figure
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=grid_plot.m>

function grid_plot(x,y,z,vorn)

% Run the tessellation
if ~exist('vorn','var')
    [~,~,vorn]=voronoisphere([x'; y'; z']);
end

% Check consistency
grumble(x,y,z,vorn);

% Plot the points
plot3(x,y,z,'k.','MarkerSize',3);
xlim([-1.1 +1.1]); ylim([-1.1 +1.1]);
zlim([-1.1 +1.1]); hold on;

% Plot the tessera
for k=1:numel(vorn)
    patch(vorn{k}(1,:),...
          vorn{k}(2,:),...
          vorn{k}(3,:),'white','FaceAlpha',1);
end

% Residual cosmetics
axis square; box on; campos([0 0 10]);
xticks([]); yticks([]); zticks([]);

end

% Consistency enforcement
function grumble(x,y,z,vorn)
if (~isnumeric(x))||(~isreal(x))||any(~isfinite(x),'all')||...
   (~isnumeric(y))||(~isreal(y))||any(~isfinite(y),'all')||...
   (~isnumeric(z))||(~isreal(z))||any(~isfinite(z),'all')
    error('all elements of x,y,z must be numeric real, and finite.');
end
if (~iscolumn(x))||(~iscolumn(y))||(~iscolumn(z))||...
   any(size(x)~=size(y),'all')||any(size(y)~=size(z),'all')
    error('x,y,z must be column vectors of the same size.');
end
if ~iscell(vorn)
    error('vorn must be a cell array of Voronoi tessera.');
end
end

% С чего начинаются линуксы?
% Со Слаки, которая труъ.
% С процесса снесения виндуса,
% И глаз, что красны поутру.
% 
% А может, они начинаются
% С болванки, что гуру принес?
% С исошек, что долго качаются,
% Щетины и длинных волос.
% 
% С чего начинаются линуксы?
% С пароля аккаунта root.
% С дурацких постов анонимусов,
% Что вечно на чанах живут.
% 
% А может, они начинаются
% С копания в гугле всегда?
% С курения документации
% И FAQ-ов: без них никуда.
% 
% С чего начинаются линуксы?
% С ядра, GPL-а и GNU,
% С консольных команд многочисленных,
% Все помнят хотя бы одну.
% 
% А может, они начинаются
% С тех Hортона синих полос?
% И клятвы которую в юности
% Ты им своим сердцем принёс.
%
% Russian Internet folklore, ca. 2005

