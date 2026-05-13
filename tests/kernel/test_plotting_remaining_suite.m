% Tests remaining Spinach plotting helper gaps under offscreen graphics. Syntax:
%
%                    result=test_plotting_remaining_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises deterministic axis, contour, cropping, molecular,
% tensor-display, ultrafast, and guarded interactive plotting helpers under
% invisible figures without relying on image comparison.
%
% ilya.kuprov@weizmann.ac.il

function result=test_plotting_remaining_suite()

% Announce the test target
fprintf('TESTING: Remaining offscreen plotting helpers\n');

% State the remaining plotting-helper target of the test
result=new_test_result('kernel/plotting_remaining_suite',...
                       'Remaining plotting helper gaps',...
                       'remaining plotting helpers must return deterministic arrays, graphics objects, or guarded validation errors.');

% Force invisible figures during the test
old_visibility=get(groot,'defaultFigureVisible');
warning_state=warning('query','MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
set(groot,'defaultFigureVisible','off');
cleaner=onCleanup(@()local_cleanup(old_visibility,warning_state));

% Build a minimal spin-system structure used by plotting routines
spin_system=local_plot_system();

% Exercise deterministic one-dimensional and transform axes
result=local_test_axes(result,spin_system);

% Exercise colour maps, contour levels, cropping, and volume zooming
result=local_test_arrays(result,spin_system);

% Exercise ordinary graphics helpers under invisible figures
result=local_test_graphics(result,spin_system);

% Exercise tensor-display graphics under invisible figures
result=local_test_tensors(result);

% Exercise non-interactive integration and guarded interactive paths
result=local_test_interactive_guards(result,spin_system);

end


function result=local_test_axes(result,spin_system)

% Check even-point Fourier axis folding
axis_even=ft_axis(1,8,4);
result=test_close(result,'ft_axis even points',axis_even,[-3 -1 1 3],1e-12,1e-12,...
                  'even Fourier axes drop the folded edge frequency');

% Check odd-point Fourier axis centring
axis_odd=ft_axis(2,10,5);
result=test_close(result,'ft_axis odd points',axis_odd,[-2 0 2 4 6],1e-12,1e-12,...
                  'odd Fourier axes are centred between the sampled edge points');

% Check sweep-to-ticks conversion
sweep_axis=sweep2ticks(10,8,5);
result=test_close(result,'sweep2ticks axis',sweep_axis,[14;12;10;8;6],1e-12,1e-12,...
                  'sweep2ticks returns descending NMR-style column ticks in Hz');

% Check FFT frequency axes with zero filling
[f_shift,f_axis,df]=fft_freq_axis(4,0.5,2);
result=test_close(result,'fft_freq_axis df',df,1/3,1e-12,1e-12,...
                  'FFT frequency resolution is sampling frequency divided by zero-filled length');
result=test_close(result,'fft_freq_axis unshifted',f_axis,(0:5)'/3,1e-12,1e-12,...
                  'unshifted FFT axes start at zero and increase by df');
result=test_close(result,'fft_freq_axis shifted',f_shift,(-3:2)'/3,1e-12,1e-12,...
                  'shifted FFT axes place negative frequencies before non-negative frequencies');

% Check IFFT time axes with symmetric zero filling
[t_shift,t_axis,dt]=ifft_time_axis(4,0.5,1);
result=test_close(result,'ifft_time_axis dt',dt,1/3,1e-12,1e-12,...
                  'IFFT time resolution is reciprocal spectral width');
result=test_close(result,'ifft_time_axis unshifted',t_axis,(0:5)'/3,1e-12,1e-12,...
                  'unshifted IFFT axes start at zero and increase by dt');
result=test_close(result,'ifft_time_axis shifted',t_shift,(-3:2)'/3,1e-12,1e-12,...
                  'shifted IFFT axes are centred around zero time');

% Check one-dimensional plotting axis in raw frequency units
parameters.spins={'1H'};
parameters.offset=100;
parameters.sweep=800;
parameters.zerofill=8;
parameters.axis_units='Hz';
[axis_hz,axis_label]=axis_1d(spin_system,parameters);
result=test_close(result,'axis_1d Hz values',axis_hz,[-300 -200 -100 0 100 200 300 400],1e-12,1e-12,...
                  'axis_1d delegates scalar sweep specifications to ft_axis');
result=test_true(result,'axis_1d Hz label',contains(axis_label,'offset frequency / Hz'),...
                 'axis_1d labels raw rotating-frame frequency axes');

% Check one-dimensional plotting axis in point units
parameters.axis_units='points';
[axis_points,axis_label]=axis_1d(spin_system,parameters);
result=test_close(result,'axis_1d point values',axis_points,1:parameters.zerofill,1e-12,1e-12,...
                  'axis_1d returns digitisation indices for point axes');
result=test_true(result,'axis_1d point label',strcmp(axis_label,'digitisation points'),...
                 'axis_1d labels point-count axes explicitly');

% Check explicit two-bound sweep specification
range_params.spins={'13C'};
range_params.sweep=[-200 300];
range_params.zerofill=6;
range_params.axis_units='Hz';
[range_axis,range_label]=axis_1d(spin_system,range_params);
result=test_close(result,'axis_1d range values',range_axis,linspace(-200,300,6),1e-12,1e-12,...
                  'axis_1d preserves explicit two-bound sweep specifications');
result=test_true(result,'axis_1d range label',contains(range_label,'offset frequency / Hz'),...
                 'axis_1d still labels explicit Hz ranges as offset frequencies');

end


function result=local_test_arrays(result,spin_system)

% Check blue-white-red colour map anchors
colour_map=bwr_cmap();
result=test_true(result,'bwr_cmap size',isequal(size(colour_map),[255 3]),...
                 'bwr_cmap returns a 255-by-3 MATLAB colour map');
result=test_close(result,'bwr_cmap blue end',colour_map(1,:),[0 0 1],1e-12,1e-12,...
                  'bwr_cmap starts at pure blue');
result=test_close(result,'bwr_cmap white centre',colour_map(128,:),[1 1 1],1e-12,1e-12,...
                  'bwr_cmap maps zero to white at the centre row');
result=test_close(result,'bwr_cmap red end',colour_map(end,:),[1 0 0],1e-12,1e-12,...
                  'bwr_cmap ends at pure red');

% Check adaptive contour spacing for mixed-sign data
[all_conts,pos_conts,neg_conts]=contspacing(10,-4,[0.1 0.3 0.2 0.4],2,'both',4);
pos_ref=(0.3-0.1)*10*linspace(0,1,4).^2+10*0.1;
neg_ref=(0.4-0.2)*(-4)*linspace(0,1,4).^2+(-4)*0.2;
result=test_close(result,'contspacing positive',pos_conts,pos_ref,1e-12,1e-12,...
                  'positive contours follow the documented non-linear spacing law');
result=test_close(result,'contspacing negative',neg_conts,neg_ref,1e-12,1e-12,...
                  'negative contours follow the documented non-linear spacing law');
result=test_close(result,'contspacing merged',all_conts,[neg_ref(end:-1:1) pos_ref],1e-12,1e-12,...
                  'merged contours list negative levels from outside in, followed by positives');

% Check positive-only contour requests
[all_pos,pos_only,neg_only]=contspacing(5,-3,[0.1 0.2 0.1 0.2],1,'positive',3);
result=test_true(result,'contspacing no negative',isempty(neg_only),...
                 'positive-only contour requests suppress negative contours');
result=test_close(result,'contspacing positive-only merge',all_pos,pos_only,1e-12,1e-12,...
                  'positive-only merged contours equal the positive contour list');

% Build deterministic two-dimensional data for cropping
spec=reshape(1:80,[8 10]);
parameters.spins={'1H','13C'};
parameters.offset=[0 0];
parameters.sweep=[1000 800];
crop_ranges={[-0.5 0.3],[-1.0 0.4]};
[crop_spec,crop_params]=crop_2d(spin_system,spec,parameters,crop_ranges);

% Compute the expected crop indices independently from the same public inputs
axis_f1_hz=linspace(-parameters.sweep(1)/2,parameters.sweep(1)/2,size(spec,1))+parameters.offset(1);
axis_f2_hz=linspace(-parameters.sweep(2)/2,parameters.sweep(2)/2,size(spec,2))+parameters.offset(2);
axis_f1_ppm=1e6*(2*pi)*axis_f1_hz/(spin(parameters.spins{1})*spin_system.inter.magnet);
axis_f2_ppm=1e6*(2*pi)*axis_f2_hz/(spin(parameters.spins{2})*spin_system.inter.magnet);
left_f1=find(axis_f1_ppm>crop_ranges{1}(1),1);
right_f1=find(axis_f1_ppm>crop_ranges{1}(2),1);
left_f2=find(axis_f2_ppm>crop_ranges{2}(1),1);
right_f2=find(axis_f2_ppm>crop_ranges{2}(2),1);
ref_spec=spec(left_f1:right_f1,left_f2:right_f2);
ref_offset=[mean(axis_f1_hz([left_f1 right_f1])) mean(axis_f2_hz([left_f2 right_f2]))];
ref_sweep=[axis_f1_hz(right_f1)-axis_f1_hz(left_f1) axis_f2_hz(right_f2)-axis_f2_hz(left_f2)];
ref_zerofill=[right_f1-left_f1 right_f2-left_f2];
result=test_close(result,'crop_2d spectrum',crop_spec,ref_spec,1e-12,1e-12,...
                  'crop_2d extracts the digital bins selected by ppm bounds');
result=test_close(result,'crop_2d offset',crop_params.offset,ref_offset,1e-12,1e-12,...
                  'crop_2d recentres offsets on the retained digital range');
result=test_close(result,'crop_2d sweep',crop_params.sweep,ref_sweep,1e-12,1e-12,...
                  'crop_2d updates sweep widths to the retained digital range');
result=test_close(result,'crop_2d zerofill',crop_params.zerofill,ref_zerofill,1e-12,1e-12,...
                  'crop_2d updates point counts to the retained digital range');

% Check three-dimensional fractional zooming
volume=reshape(1:1000,[10 10 10]);
ext=[0 9 10 19 -5 4];
zoom_ranges=[0.2 0.6 0.3 0.8 0.1 0.5];
[zoom_volume,zoom_ext]=zoom_3d(volume,ext,zoom_ranges);
result=test_close(result,'zoom_3d volume',zoom_volume,volume(2:6,3:8,1:5),1e-12,1e-12,...
                  'zoom_3d extracts the subcube implied by fractional bounds');
result=test_close(result,'zoom_3d extents',zoom_ext,[1 5 12 17 -5 -1],1e-12,1e-12,...
                  'zoom_3d maps retained indices back to physical extents');

end


function result=local_test_graphics(result,spin_system)

% Exercise molecular stick plotting with supplied connectivity
fig=figure('Visible','off');
xyz=[0 0 0;1 0 0;1 1 0];
conmatrix=logical([0 1 0;0 0 1;0 0 0]);
molplot(xyz,conmatrix);
line_obj=findobj(fig,'Type','line');
x_data=get(line_obj(1),'XData');
result=test_true(result,'molplot line count',isscalar(line_obj),...
                 'molplot renders supplied connectivity as one line object with NaN-separated bonds');
result=test_true(result,'molplot bond separators',nnz(isnan(x_data))==2,...
                 'molplot separates two supplied bonds by NaN coordinates');
close(fig);

% Exercise cylindrical grid drawing under invisible graphics
fig=figure('Visible','off');
cylgrid(-1,2,3);
line_obj=findobj(fig,'Type','line');
text_obj=findobj(fig,'Type','text');
axes_obj=findobj(fig,'Type','axes');
result=test_true(result,'cylgrid line objects',numel(line_obj)>=57,...
                 'cylgrid draws spokes, rings, and vertical tick marks');
result=test_true(result,'cylgrid text objects',numel(text_obj)>=19,...
                 'cylgrid draws azimuth and vertical tick labels');
result=test_true(result,'cylgrid hidden axes',strcmp(get(axes_obj(1),'Visible'),'off'),...
                 'cylgrid hides default Cartesian axes after drawing the grid');
close(fig);

% Exercise ultrafast 2D plotting with consistent compact dimensions
fig=figure('Visible','off');
uf_params=local_uf_params();
uf_rows=local_uf_rows(uf_params);
[row_grid,col_grid]=ndgrid(linspace(-1,1,uf_rows),linspace(-1,1,uf_params.nloops));
uf_spectrum=exp(-row_grid.^2-col_grid.^2)-0.2*exp(-6*((row_grid-0.3).^2+(col_grid+0.2).^2));
plot_uf(spin_system,uf_spectrum,uf_params);
contour_obj=findobj(fig,'Type','contour');
axes_obj=findobj(fig,'Type','axes');
result=test_true(result,'plot_uf contour object',isscalar(contour_obj),...
                 'plot_uf creates one contour object for compact ultrafast spectra');
result=test_true(result,'plot_uf reversed axes',strcmp(get(axes_obj(1),'XDir'),'reverse')&&...
                 strcmp(get(axes_obj(1),'YDir'),'reverse'),...
                 'plot_uf applies NMR-style reversed axes in both dimensions');
close(fig);

end


function result=local_test_tensors(result)

% Build minimal tensor-display input structures
props=local_tensor_props();
options.style='ellipsoids';
options.kill_iso=false;
options.numbers=false;
options.symbols=false;
conmatrix=false(2);

% Exercise chemical shielding tensor display
fig=figure('Visible','off');
cst_display(props,1,0.05,conmatrix,options);
surf_obj=findobj(fig,'Type','surface');
line_obj=findobj(fig,'Type','line');
result=test_true(result,'cst_display surface',numel(surf_obj)>=1,...
                 'cst_display draws at least one tensor ellipsoid surface');
result=test_true(result,'cst_display eigenvectors',numel(line_obj)>=3,...
                 'cst_display draws tensor eigensystem vectors as line objects');
close(fig);

% Exercise electric field gradient tensor display
fig=figure('Visible','off');
efg_display(props,1,0.05,conmatrix,options);
surf_obj=findobj(fig,'Type','surface');
line_obj=findobj(fig,'Type','line');
result=test_true(result,'efg_display surface',numel(surf_obj)>=1,...
                 'efg_display draws at least one tensor ellipsoid surface');
result=test_true(result,'efg_display eigenvectors',numel(line_obj)>=3,...
                 'efg_display draws tensor eigensystem vectors as line objects');
close(fig);

% Exercise hyperfine coupling tensor display
fig=figure('Visible','off');
hfc_display(props,1,0.05,conmatrix,options);
surf_obj=findobj(fig,'Type','surface');
line_obj=findobj(fig,'Type','line');
result=test_true(result,'hfc_display surface',numel(surf_obj)>=1,...
                 'hfc_display draws at least one tensor ellipsoid surface');
result=test_true(result,'hfc_display eigenvectors',numel(line_obj)>=3,...
                 'hfc_display draws tensor eigensystem vectors as line objects');
close(fig);

end


function result=local_test_interactive_guards(result,spin_system)

% Build deterministic data for non-interactive 2D integration
parameters.spins={'1H','13C'};
parameters.offset=[0 0];
parameters.sweep=[1000 800];
parameters.zerofill=[8 10];
parameters.axis_units='Hz';
[row_grid,col_grid]=ndgrid(linspace(-1,1,8),linspace(-1,1,10));
spectrum=exp(-row_grid.^2-col_grid.^2)-0.3*exp(-8*((row_grid-0.4).^2+(col_grid+0.3).^2));

% Save a temporary integration-range file to avoid mouse input
range_file=[tempname '.mat'];
ranges={[-100;100],[-80;80]};
save(range_file,'ranges');
file_cleaner=onCleanup(@()local_delete_file(range_file));

% Exercise int_2d through the file-driven path
fig=figure('Visible','off');
int_2d(spin_system,spectrum,parameters,4,[0.1 0.8 0.1 0.8],1,32,2,'both',range_file);
contour_obj=findobj(fig,'Type','contour');
result=test_true(result,'int_2d contour object',isscalar(contour_obj),...
                 'int_2d reuses plot_2d while integrating ranges from a file without mouse input');
result=test_true(result,'int_2d range file retained',exist(range_file,'file')==2,...
                 'int_2d consumes an existing range file without replacing it interactively');
close(fig);

% Guard slice_2d by checking its validation path before mouse input
bad_call=@()slice_2d(spin_system,spectrum,parameters,0,[0.1 0.8 0.1 0.8],1,32,2,'both');
result=local_expect_error(result,'slice_2d invalid ncont',bad_call,'ncont parameter must be a positive integer',...
                          'slice_2d rejects invalid contour counts before entering its mouse-driven loop');
result.messages{end+1}='SKIP: slice_2d full extraction loop requires live ginput interaction and is intentionally not driven in offscreen automation.';

% Guard write_movie cheaply through input validation
movie_bad_call=@()write_movie(17);
result=local_expect_error(result,'write_movie invalid filename',movie_bad_call,'file_name must be a character string',...
                          'write_movie rejects non-character filenames before opening VideoWriter');

% Optionally exercise the slow movie writer only when explicitly requested
if strcmp(getenv('SPINACH_RUN_SLOW_PLOTTING'),'1')
    result=local_test_write_movie(result);
else
    result.messages{end+1}='SKIP: write_movie full MP4 generation is slow, codec-dependent, and gated by SPINACH_RUN_SLOW_PLOTTING=1.';
end

end


function result=local_test_write_movie(result)

% Prepare a tiny three-dimensional figure for movie capture
fig=figure('Visible','off');
plot3([0 1],[0 1],[0 1],'k-');
axis vis3d;
movie_file=[tempname '.mp4'];
file_cleaner=onCleanup(@()local_delete_file(movie_file));

% Run the movie writer and check that it produced data
write_movie(movie_file);
movie_info=dir(movie_file);
result=test_true(result,'write_movie file created',exist(movie_file,'file')==2&&movie_info.bytes>0,...
                 'write_movie writes a non-empty MPEG-4 file when the graphics backend and codec are available');
close(fig);

end


function result=local_expect_error(result,label,call_handle,message_part,why)

% Run the supplied call and require a matching failure
try
    call_handle();
    got_error=false;
    err_msg='';
catch err
    got_error=true;
    err_msg=err.message;
end

% Check that the failure was the expected validation error
result=test_true(result,label,got_error&&contains(err_msg,message_part),why);

end


function spin_system=local_plot_system()

% Provide the minimal fields used by plotting utilities
spin_system.comp.isotopes={'1H','13C'};
spin_system.inter.magnet=14.1;
spin_system.sys.disable={};
spin_system.sys.output='hush';
spin_system.tols.freeg=2.00231930436256;

end


function props=local_tensor_props()

% Provide compact geometry and tensor data for display helpers
props.std_geom=[0 0 0;1 0 0];
props.symbols={'C','H'};
props.cst={diag([1 2 3]),diag([-1 0.5 2])};
props.efg={diag([1 -2 1]),diag([2 -1 -1])};
props.nqi={[],[]};
props.hfc.full.matrix={diag([1 2 3]),diag([-1 2 4])};

end


function parameters=local_uf_params()

% Provide compact ultrafast plotting parameters with integer k-space rows
parameters.spins={'1H','13C'};
parameters.dims=0.002;
parameters.deltat=1e-4;
parameters.npoints=10;
parameters.Ga=0.1;
parameters.offset=[0 0];
parameters.axis_units='Hz';
parameters.offset_uf_cov=0;
parameters.nloops=8;
parameters.Te=1e-3;

end


function rows=local_uf_rows(parameters)

% Compute the required ultrafast dimension size from plot_uf's public formula
Ta=parameters.deltat*parameters.npoints;
gamma=spin(parameters.spins{1});
k_max=gamma*parameters.Ga*Ta/(2*pi);
rows=round(parameters.dims*k_max);

end


function local_delete_file(file_name)

% Delete a temporary file if it survived the test body
if exist(file_name,'file')
    delete(file_name);
end

end


function local_cleanup(old_visibility,warning_state)

% Restore figure and warning state after success or failure
close all force;
set(groot,'defaultFigureVisible',old_visibility);
warning(warning_state.state,warning_state.identifier);

end


