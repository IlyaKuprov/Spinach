% Tests offscreen execution of Spinach plotting helpers. Syntax:
%
%                    result=test_plotting_helpers_offscreen()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises plotting helpers under invisible figures, and checks
% graphics object creation, axis sizes, returned data arrays, and figure
% helper side effects without relying on image comparison.
%
% ilya.kuprov@weizmann.ac.il

function result=test_plotting_helpers_offscreen()

% Announce the test target
fprintf('TESTING: Offscreen plotting helpers\n');

% State the plotting-helper target of the test
result=new_test_result('kernel/plotting_helpers_offscreen',...
                       'Offscreen plotting helpers',...
                       'plotting helpers must create deterministic graphics objects under invisible offscreen figures.');

% Force invisible figures during the test
old_visibility=get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
cleaner=onCleanup(@()local_cleanup(old_visibility));

% Build a minimal spin-system structure used only by plotting routines
spin_system=local_plot_system();

% Exercise house-style figure helpers
result=local_test_house_style(result);

% Exercise one-dimensional spectral plotting
result=local_test_plot_1d(result,spin_system);

% Exercise two-dimensional contour and stack plotting
result=local_test_plot_2d(result,spin_system);

% Exercise three-dimensional spectrum plotting
result=local_test_plot_3d(result,spin_system);

% Exercise Bloch-sphere, MRI, and volume plotting utilities
result=local_test_misc_plots(result);

end


function result=local_test_house_style(result)

% Create an invisible figure and check default handling
fig=kfigure('Visible','off');
result=test_true(result,'kfigure invisible handle',ishandle(fig)&&strcmp(get(fig,'Visible'),'off'),...
                 'kfigure returns a valid invisible figure handle when requested');

% Scale the figure and check deterministic dimensions
scale_figure([1.25 0.75]);
fig_pos=get(fig,'Position');
def_pos=get(groot,'defaultFigurePosition');
result=test_close(result,'scale_figure width',fig_pos(3),1.25*def_pos(3),1e-12,1e-12,...
                  'scale_figure scales width from the root default figure size');
result=test_close(result,'scale_figure height',fig_pos(4),0.75*def_pos(4),1e-12,1e-12,...
                  'scale_figure scales height from the root default figure size');

% Exercise subplot and labelling helpers on a 3D line
ksubplot(1,2,1);
plot3([0 1],[0 1],[0 1],'k-');
kgrid(); ktitle('helper title');
kxlabel('X axis'); kylabel('Y axis'); kzlabel('Z axis');
leg_obj=klegend({'trace'});

% Exercise super-title and colour-bar helpers
ksubplot(1,2,2);
imagesc(magic(3));
kcolourbar('intensity');
ksgtitle('helper grid');

% Check that helper calls populated graphics metadata
axes_obj=findall(fig,'Type','axes');
colour_obj=findall(fig,'Type','colorbar');
result=test_true(result,'ksubplot axes count',numel(axes_obj)>=2,...
                 'ksubplot creates ordinary axes objects');
result=test_true(result,'klegend object',ishandle(leg_obj),...
                 'klegend returns a valid legend object');
result=test_true(result,'kcolourbar object',isscalar(colour_obj),...
                 'kcolourbar creates one colour-bar object');
close(fig);

end


function result=local_test_plot_1d(result,spin_system)

% Define a deterministic one-dimensional spectrum
parameters.spins={'1H'};
parameters.offset=0;
parameters.sweep=1000;
parameters.zerofill=32;
parameters.axis_units='Hz';
parameters.invert_axis=0;
x_axis=linspace(-2,2,parameters.zerofill)';
spectrum=exp(-x_axis.^2);

% Plot the spectrum under an invisible figure
fig=kfigure('Visible','off');
plot_1d(spin_system,spectrum,parameters,'k-');
line_obj=findobj(fig,'Type','line');
axes_obj=findobj(fig,'Type','axes');
y_data=get(line_obj(1),'YData')';

% Check that the plotted curve matches the input data
result=test_true(result,'plot_1d line count',isscalar(line_obj),...
                 'real one-dimensional spectra produce one line object');
result=test_close(result,'plot_1d y data',y_data,spectrum,1e-12,1e-12,...
                  'plot_1d forwards the real spectrum to Matlab line data');
result=test_true(result,'plot_1d axis direction',strcmp(get(axes_obj(1),'XDir'),'normal'),...
                 'plot_1d respects parameters.invert_axis=0');
close(fig);

% Plot a complex spectrum to exercise recursive real and imaginary handling
fig=kfigure('Visible','off');
plot_1d(spin_system,spectrum+1i*(2*spectrum),parameters,'k-');
line_obj=findobj(fig,'Type','line');
legend_obj=findobj(fig,'Type','legend');
result=test_true(result,'plot_1d complex line count',numel(line_obj)==2,...
                 'complex one-dimensional spectra produce real and imaginary line objects');
result=test_true(result,'plot_1d complex legend',isscalar(legend_obj),...
                 'complex one-dimensional spectra create a real/imaginary legend');
close(fig);

end


function result=local_test_plot_2d(result,spin_system)

% Define a deterministic two-dimensional spectrum
parameters.spins={'1H','13C'};
parameters.offset=[0 0];
parameters.sweep=[1000 800];
parameters.zerofill=[12 16];
parameters.axis_units='Hz';
spectrum=local_spectrum_2d(parameters.zerofill(1),parameters.zerofill(2));

% Contour-plot the two-dimensional spectrum
fig=kfigure('Visible','off');
[axis_f1,axis_f2,plot_spectrum]=plot_2d(spin_system,spectrum,parameters,...
                                        4,[0.10 0.80 0.10 0.80],1,32,2,'both');
contour_obj=findobj(fig,'Type','contour');
colour_obj=findobj(fig,'Type','colorbar');
axes_obj=findobj(fig,'Type','axes');

% Check returned axes, transposed data, and graphics objects
result=test_close(result,'plot_2d returned spectrum',plot_spectrum,transpose(spectrum),1e-12,1e-12,...
                  'plot_2d returns the transposed spectrum used for contour plotting');
result=test_true(result,'plot_2d axis sizes',numel(axis_f1)==size(spectrum,2)&&numel(axis_f2)==size(spectrum,1),...
                 'plot_2d returns frequency axes matching the plotted matrix dimensions');
result=test_true(result,'plot_2d contour object',isscalar(contour_obj),...
                 'plot_2d creates a contour object');
result=test_true(result,'plot_2d colour bar',isscalar(colour_obj),...
                 'plot_2d creates one colour bar unless disabled');
result=test_true(result,'plot_2d reversed axes',strcmp(get(axes_obj(1),'XDir'),'reverse')&&strcmp(get(axes_obj(1),'YDir'),'reverse'),...
                 'plot_2d uses NMR-style reversed axes');
close(fig);

% Stack-plot a square spectrum along the indirect dimension
stack_params=parameters;
stack_params.zerofill=[12 12];
stack_spectrum=local_spectrum_2d(stack_params.zerofill(1),stack_params.zerofill(2));
fig=kfigure('Visible','off');
stack_2d(spin_system,stack_spectrum,stack_params,1,@(slice)norm(slice,2));
patch_obj=findobj(fig,'Type','patch');
patch_x=get(patch_obj(1),'XData');
result=test_true(result,'stack_2d patch count',numel(patch_obj)==size(stack_spectrum,2),...
                 'stack_2d creates one patch line per selected spectral slice');
result=test_true(result,'stack_2d patch length',numel(patch_x)==size(stack_spectrum,1)+1,...
                 'stack_2d closes each patch line with a trailing NaN point');
close(fig);

end


function result=local_test_plot_3d(result,spin_system)

% Define a compact deterministic three-dimensional spectrum
parameters.spins={'1H','13C','15N'};
parameters.offset=[0 0 0];
parameters.sweep=[900 700 500];
parameters.npoints=[5 5 5];
parameters.zerofill=[5 5 5];
parameters.axis_units='Hz';
[x_grid,y_grid,z_grid]=ndgrid(linspace(-1,1,5),linspace(-1,1,5),linspace(-1,1,5));
spectrum=exp(-(x_grid.^2+y_grid.^2+z_grid.^2));
spectrum=spectrum-mean(spectrum(:));

% Plot isosurfaces and the three projections
fig=kfigure('Visible','off');
plot_3d(spin_system,spectrum,parameters,2,[0.20 0.80 0.20 0.80],1,'both');
patch_obj=findobj(fig,'Type','patch');
axes_obj=findobj(fig,'Type','axes');

% Check that volume and projection graphics were created
result=test_true(result,'plot_3d patch objects',numel(patch_obj)>=2,...
                 'plot_3d creates patch objects for isosurfaces and projection contours');
result=test_true(result,'plot_3d axes count',numel(axes_obj)>=4,...
                 'plot_3d creates the volume view and three projection axes');
close(fig);

end


function result=local_test_misc_plots(result)

% Exercise Bloch-sphere trajectory plotting
x_traj=cos(linspace(0,pi,16));
y_traj=sin(linspace(0,pi,16));
z_traj=linspace(-1,1,16);
bloch_sph_plot(x_traj,y_traj,z_traj,'k-');
fig=gcf();
surf_obj=findobj(fig,'Type','surface');
line_obj=findobj(fig,'Type','line');
axis_obj=findobj(fig,'Type','line','Tag','InstantaneousAxis');
result=test_true(result,'bloch_sph_plot surface',numel(surf_obj)>=1,...
                 'bloch_sph_plot creates a reference sphere surface');
result=test_true(result,'bloch_sph_plot trajectory',numel(line_obj)>=1,...
                 'bloch_sph_plot creates a trajectory line');
result=test_true(result,'bloch_sph_plot axis curve',numel(axis_obj)==1,...
                 'bloch_sph_plot creates an instantaneous rotation-axis curve');
axis_xyz=[get(axis_obj,'XData'); get(axis_obj,'YData'); get(axis_obj,'ZData')];
axis_rad=sqrt(sum(axis_xyz.^2,1));
result=test_true(result,'bloch_sph_plot axis normalisation',...
                 all(axis_rad(isfinite(axis_rad))<=1+1e-12),...
                 'bloch_sph_plot keeps the instantaneous axis inside the unit sphere');
close(fig);

% Check row-wise static plotting of multiple trajectories
multi_time=linspace(0,pi/2,5);
x_multi=[cos(multi_time); cos(multi_time)];
y_multi=[sin(multi_time); -sin(multi_time)];
z_multi=zeros(size(x_multi));
bloch_sph_plot(x_multi,y_multi,z_multi,'k-');
fig=gcf();
axis_obj=findobj(fig,'Type','line','Tag','InstantaneousAxis');
axis_zdata=[];
if isscalar(axis_obj)
    axis_zdata=get(axis_obj,'ZData');
end
npoints=numel(multi_time);
result=test_true(result,'bloch_sph_plot static multirow axis',...
                 isscalar(axis_obj)&&numel(axis_zdata)==2*(npoints+1)&&...
                 all(axis_zdata(1:npoints)>0)&&isnan(axis_zdata(npoints+1))&&...
                 all(axis_zdata((npoints+2):(2*npoints+1))<0)&&isnan(axis_zdata(end)),...
                 'bloch_sph_plot traces each row trajectory separately in static instantaneous-axis plotting');
close(fig);

% Optionally exercise the Bloch-sphere movie path
if strcmp(getenv('SPINACH_RUN_SLOW_PLOTTING'),'1')
    movie_file=[tempname '.avi'];
    x_movie=cos(linspace(0,pi/2,5));
    y_movie=sin(linspace(0,pi/2,5));
    z_movie=zeros(size(x_movie));
    bloch_sph_plot(x_movie,y_movie,z_movie,'k-',movie_file);
    fig=gcf();
    axis_obj=findobj(fig,'Type','line','Tag','InstantaneousAxis');
    axis_xdata=[];
    if isscalar(axis_obj)
        axis_xdata=get(axis_obj,'XData');
    end
    movie_exists=(exist(movie_file,'file')==2);
    movie_bytes=0;
    if movie_exists
        movie_info=dir(movie_file);
        movie_bytes=movie_info.bytes;
    end
    result=test_true(result,'bloch_sph_plot movie file',...
                     movie_exists&&movie_bytes>0,...
                     'bloch_sph_plot writes a non-empty video file in movie mode');
    result=test_true(result,'bloch_sph_plot movie axis stick',...
                     isscalar(axis_obj)&&numel(axis_xdata)==3&&...
                     axis_xdata(1)==0&&isnan(axis_xdata(end)),...
                     'bloch_sph_plot leaves a moving instantaneous-axis stick, not an axis-tip curve, in movie mode');
    frame_count=0;
    if movie_exists
        movie_reader=VideoReader(movie_file);
        while hasFrame(movie_reader)
            readFrame(movie_reader);
            frame_count=frame_count+1;
        end
    end
    result=test_true(result,'bloch_sph_plot movie frame count',frame_count==numel(x_movie),...
                     'bloch_sph_plot writes one movie frame per input trajectory point');
    close(fig);
    if exist(movie_file,'file')==2
        delete(movie_file);
    end
else
    result.messages{end+1}='SKIP: bloch_sph_plot video generation is codec-dependent and gated by SPINACH_RUN_SLOW_PLOTTING=1.';
end

% Exercise MRI image and k-space plotting branches
fig=kfigure('Visible','off');
parameters.spins={'1H'};
parameters.pe_grad_amp=0.01;
parameters.pe_grad_dur=1e-3;
parameters.ro_grad_amp=0.02;
parameters.ro_grad_dur=1e-3;
parameters.dims=[0.02 0.03];
mri_image=reshape(1:16,[4 4]);
subplot(1,3,1); mri_2d_plot(mri_image,parameters,'image');
subplot(1,3,2); mri_2d_plot(mri_image,parameters,'phantom');
subplot(1,3,3); mri_2d_plot(mri_image+1i*mri_image,parameters,'k-space');
image_obj=findobj(fig,'Type','image');
result=test_true(result,'mri_2d_plot image objects',numel(image_obj)==3,...
                 'mri_2d_plot creates one image object for each supported plotting branch');
result=test_close(result,'mri_2d_plot k-space real data',get(image_obj(1),'CData'),real(mri_image+1i*mri_image),1e-12,1e-12,...
                  'mri_2d_plot displays the real part of complex k-space data');
close(fig);

% Exercise volumetric plotting on a compact signed cube
fig=kfigure('Visible','off');
data_cube=zeros(4,4,4);
data_cube(2,2,2)=1;
data_cube(3,3,3)=-0.5;
volplot(data_cube,[-1 1 -1 1 -1 1],[1 1]);
surf_obj=findobj(fig,'Type','surface');
result=test_true(result,'volplot surfaces',numel(surf_obj)>=1,...
                 'volplot creates semitransparent surfaces for non-zero volume slices');
close(fig);

end


function spin_system=local_plot_system()

% Provide the minimal fields used by plotting utilities
spin_system.comp.isotopes={'1H','13C','15N'};
spin_system.inter.magnet=14.1;
spin_system.sys.disable={};
spin_system.sys.output='hush';

end


function spectrum=local_spectrum_2d(nrows,ncols)

% Build smooth positive and negative peaks without random numbers
[row_grid,col_grid]=ndgrid(linspace(-1,1,nrows),linspace(-1,1,ncols));
pos_peak=exp(-8*((row_grid-0.30).^2+(col_grid+0.20).^2));
neg_peak=0.75*exp(-10*((row_grid+0.25).^2+(col_grid-0.35).^2));
spectrum=pos_peak-neg_peak;

end


function local_cleanup(old_visibility)

% Restore figure state after success or failure
close all force;
set(groot,'defaultFigureVisible',old_visibility);

end
