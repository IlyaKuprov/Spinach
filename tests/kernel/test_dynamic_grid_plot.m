% Tests grid_plot() under offscreen graphics. Syntax:
%
%                    result=test_dynamic_grid_plot()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test draws a tetrahedral spherical Voronoi tessellation, checks the
% patch count and numeric colour mapping, and verifies optional centre dots.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_grid_plot()

% Announce the test target
fprintf('TESTING: Spherical grid plotting\n');

% State the grid_plot target of the test
result=new_test_result('kernel/dynamic_grid_plot',...
                       'Offscreen spherical grid plotting',...
                       'grid_plot() must draw one patch per Voronoi cell and honour centre-dot options.');

% Force invisible figures during the test
old_visibility=get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
cleaner=onCleanup(@()set(groot,'defaultFigureVisible',old_visibility));

% Build a regular tetrahedral grid on the unit sphere
xyz=[1 1 -1 -1;1 -1 1 -1;1 -1 -1 1];
xyz=xyz./sqrt(sum(xyz.^2,1));
x=xyz(1,:).';
y=xyz(2,:).';
z=xyz(3,:).';

% Compute the spherical Voronoi tessellation once
[~,~,vorn]=voronoisphere(xyz,pi/3);

% Draw supplied tessera with numeric colours and no centre dots
fig=figure('Visible','off');
options.dots=false;
colours=(1:4).';
grid_plot(x,y,z,vorn,colours,options);
patch_obj=findobj(fig,'Type','patch');
line_obj=findobj(fig,'Type','line');
colour_data=local_patch_colours(patch_obj);
result=test_close(result,'grid_plot supplied patch count',numel(patch_obj),numel(vorn),0,0,...
                  'a tetrahedral Voronoi grid must render one patch per tessellation cell');
result=test_close(result,'grid_plot numeric colours',sort(colour_data),colours,1e-15,1e-15,...
                  'numeric colour data must be passed through to patch objects cell by cell');
result=test_true(result,'grid_plot dots disabled',isempty(line_obj),...
                 'setting options.dots=false must suppress centre marker plotting');
close(fig);

% Draw with internally generated tessera and default centre dots
fig=figure('Visible','off');
grid_plot(x,y,z);
patch_obj=findobj(fig,'Type','patch');
line_obj=findobj(fig,'Type','line');
axes_obj=findobj(fig,'Type','axes');
result=test_close(result,'grid_plot generated patch count',numel(patch_obj),numel(vorn),0,0,...
                  'when tessera are omitted, grid_plot() must generate the same cell count');
result=test_true(result,'grid_plot default dots',isscalar(line_obj),...
                 'the default options must draw one line object containing centre dots');
result=test_true(result,'grid_plot square axes',strcmp(get(axes_obj(1),'PlotBoxAspectRatioMode'),'manual'),...
                 'grid_plot() must leave the axes in square plotting mode');
close(fig);

end


function colour_data=local_patch_colours(patch_obj)

% Preallocate colour data
colour_data=zeros(numel(patch_obj),1);

% Pull scalar colour data out of each patch
for n=1:numel(patch_obj)
    cdata=get(patch_obj(n),'CData');
    colour_data(n)=cdata(1);
end

end


