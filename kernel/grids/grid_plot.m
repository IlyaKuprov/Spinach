% Spherical quadrature grid plotter. Takes a cloud of points
% on a sphere and plots its Voronoi tessellation. Syntax:
%
%               grid_plot(x,y,z,vorn,c,options)
%
% Parameters:
%
%   x,y,z  - column vectors containing Cartesian 
%            coordinates of grid points
%
%   c      - values to be mapped into the colour
%            of each tessellation face, white if
%            this input is left empty
%
%   vorn   - Voronoi tessellation; if this is not 
%            provided, it will be computed
%
%   options.dots - the default (true) puts black
%                  dots at centres of tessellati-
%                  on faces
%
% Outputs:
%
%   this function plots a figure
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=grid_plot.m>

function grid_plot(x,y,z,vorn,c,options)

% Run Voronoi tessellation on a sphere
if (~exist('vorn','var'))||isempty(vorn)
    [~,~,vorn]=voronoisphere([x'; y'; z']);
end

% Default tessellation face colour is white
if (~exist('c','var'))||isempty(c), c='white'; end

% Default is to draw centre dots
if (~exist('options','var'))||...
   (~isfield(options,'dots'))
    options.dots=true;
end

% Check consistency
grumble(x,y,z,vorn);

% Plot the points
if options.dots
    plot3(x,y,z,'k.','MarkerSize',3);
end
xlim([-1.1 +1.1]); ylim([-1.1 +1.1]);
zlim([-1.1 +1.1]); hold on;

% Plot the tessera
for k=1:numel(vorn)
    if ischar(c)

        % Colour by the name
        patch(vorn{k}(1,:),vorn{k}(2,:),...
              vorn{k}(3,:),c,'FaceAlpha',1);

    else

        % Colour by the numbers
        patch(vorn{k}(1,:),vorn{k}(2,:),...
              vorn{k}(3,:),c(k),'FaceAlpha',1);

    end
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

