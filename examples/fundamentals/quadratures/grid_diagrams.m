% Spherical grid diagrams for IK's book.
%
% ilya.kuprov@weizmann.ac.il

function grid_diagrams()

% Plotting logistics
kfigure(); scale_figure([1.5 1.3]);
tiledlayout(2,3,'TileSpacing','Compact',...
                'Padding','Compact');
            
% Number sequence grid example - Fibonacci
nexttile; grid_fibon('fib',321);  text(-1.0,0.95,'SEQ');
text(1.0,-0.95,'643','HorizontalAlignment','right',...
                     'Color',[0.75 0.75 0.75]);

% Polyhedron subdivision grid - icosahedral
load('../../../kernel/grids/icos_2ang_642pts.mat','betas','gammas');
nexttile; grid_plot(sin(betas).*cos(gammas),...
                    sin(betas).*sin(gammas),...
                    cos(betas)); 
text(-1.0,0.95,'ICO');
text(1.0,-0.95,'642','HorizontalAlignment','right',...
                     'Color',[0.75 0.75 0.75]);         

% Polyhedron subdivision grid - octahedral
nexttile; grid_trian('stoll',13); text(-1.0,0.95,'OCT');
text(1.0,-0.95,'678','HorizontalAlignment','right',...
                     'Color',[0.75 0.75 0.75]);
                 
% Optimisation grid - repulsion
[~,betas,gammas]=repulsion(620,3,200);
nexttile; grid_plot(sin(betas).*cos(gammas),...
                    sin(betas).*sin(gammas),...
                    cos(betas)); 
text(-1.0,0.95,'OPT');
text(1.0,-0.95,'620','HorizontalAlignment','right',...
                     'Color',[0.75 0.75 0.75]);

% Natual world inspiration grid - Igloo
nexttile; grid_igloo(23); text(-1.0,0.95,'NAT');
text(1.0,-0.95,'616','HorizontalAlignment','right',...
                     'Color',[0.75 0.75 0.75]);
                 
% "You are all wankers" - Vyacheslav Lebedev
load('../../../kernel/grids/leb_2ang_rank_41.mat','betas','gammas');
nexttile; grid_plot(sin(betas).*cos(gammas),...
                    sin(betas).*sin(gammas),...
                    cos(betas)); 
text(-1.0,0.95,'LEB');
text(1.0,-0.95,'590','HorizontalAlignment','right',...
                     'Color',[0.75 0.75 0.75]);

end

