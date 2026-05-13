% Tests compact dynamic example-stage execution with plotting. Syntax:
%
%                    result=test_dynamic_examples_smoke()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test runs short liquid-state NMR calculations adapted from plotting
% examples, processes the deterministic signals, and verifies the plotted
% graphics objects under invisible offscreen figures.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_examples_smoke()

% State the dynamic example-stage target of the test
result=new_test_result('examples/dynamic_examples_smoke',...
                       'Dynamic plotting example smoke paths',...
                       'compact example-stage calculations must run, process, and plot deterministic spectra.');

% Force invisible figures during the test
old_visibility=get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
cleaner=onCleanup(@()local_cleanup(old_visibility));

% Run and plot a compact one-dimensional acquisition path
result=local_test_acquire_1d(result);

% Run and plot a compact two-dimensional CT-COSY path
result=local_test_ct_cosy_2d(result);

end


function result=local_test_acquire_1d(result)

% Build a zero-offset one-spin Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Set up a compact free-induction acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=1000;
parameters.npoints=8;
parameters.zerofill=8;
parameters.axis_units='Hz';
parameters.invert_axis=0;

% Run the production liquid-state acquisition pathway
fid=liquid(spin_system,@acquire,parameters,'nmr');
spectrum=fftshift(fft(fid,parameters.zerofill));

% The zero-offset signal is constant and Fourier-transforms to DC only
fid_ref=0.5*ones(parameters.npoints,1);
spectrum_ref=zeros(parameters.zerofill,1);
spectrum_ref(parameters.zerofill/2+1)=parameters.zerofill*fid_ref(1);
result=test_close(result,'one-spin FID',fid,fid_ref,1e-12,1e-12,...
                  'zero-offset transverse magnetisation remains constant during acquisition');
result=test_close(result,'one-spin spectrum',spectrum,spectrum_ref,1e-12,1e-12,...
                  'the Fourier transform of a constant FID has only the centred DC point');

% Plot the deterministic one-dimensional spectrum
fig=kfigure('Visible','off');
plot_1d(spin_system,real(spectrum),parameters,'k-');
line_obj=findobj(fig,'Type','line');
y_data=get(line_obj(1),'YData')';
result=test_close(result,'one-spin plotted spectrum',y_data,real(spectrum),1e-12,1e-12,...
                  'plot_1d receives the processed one-spin spectrum from the dynamic pathway');
close(fig);

end


function result=local_test_ct_cosy_2d(result)

% Build the two-spin system used in the CT-COSY plotting example
sys.isotopes={'1H','1H'};
sys.magnet=5.9;
inter.zeeman.scalar={1.00 3.00};
inter.coupling.scalar{1,2}=7.0;
inter.coupling.scalar{2,2}=0;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Use compact point counts while preserving the production example stages
parameters.offset=500;
parameters.sweep=[2000 2000];
parameters.npoints=[8 8];
parameters.zerofill=[16 16];
parameters.spins={'1H'};
parameters.axis_units='Hz';

% Run the production CT-COSY pulse sequence through liquid()
fid=liquid(spin_system,@ct_cosy,parameters,'nmr');
spectrum=fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)));
abs_spectrum=abs(spectrum);

% Check deterministic sizes and numerical invariants from the compact run
result=test_true(result,'CT-COSY FID size',isequal(size(fid),parameters.npoints),...
                 'the compact CT-COSY pathway returns the requested time-domain matrix');
result=test_true(result,'CT-COSY spectrum size',isequal(size(spectrum),parameters.zerofill),...
                 'the compact CT-COSY pathway returns the requested zero-filled spectrum');
result=test_close(result,'CT-COSY FID norm',norm(fid(:),2),5.6290947235800708e+00,1e-10,1e-10,...
                  'the compact CT-COSY time-domain matrix is deterministic');
result=test_close(result,'CT-COSY spectrum norm',norm(spectrum(:),2),9.0065515577281133e+01,1e-10,1e-10,...
                  'the compact CT-COSY spectrum norm is deterministic');
result=test_close(result,'CT-COSY spectrum sum',sum(abs_spectrum(:)),7.1504247573723626e+02,1e-9,1e-10,...
                  'the compact CT-COSY absolute spectrum sum is deterministic');
result=test_close(result,'CT-COSY spectrum maximum',max(abs_spectrum(:)),3.0877202568062110e+01,1e-10,1e-10,...
                  'the compact CT-COSY peak amplitude is deterministic');

% Plot the processed two-dimensional spectrum
fig=kfigure('Visible','off');
scale_figure([1.5 2.0]);
[axis_f1,axis_f2,plot_spectrum]=plot_2d(spin_system,abs_spectrum,parameters,...
                                        6,[0.10 0.50 0.10 0.50],2,64,4,'positive');
contour_obj=findobj(fig,'Type','contour');
result=test_close(result,'CT-COSY plotted spectrum',plot_spectrum,transpose(abs_spectrum),1e-12,1e-12,...
                  'plot_2d receives the processed compact CT-COSY spectrum');
result=test_true(result,'CT-COSY plot axes',numel(axis_f1)==parameters.zerofill(2)&&numel(axis_f2)==parameters.zerofill(1),...
                 'plot_2d returns axes matching the compact CT-COSY zero-filled dimensions');
result=test_true(result,'CT-COSY contour object',isscalar(contour_obj),...
                 'plot_2d creates one contour object for the compact CT-COSY spectrum');
close(fig);

end


function local_cleanup(old_visibility)

% Restore figure state after success or failure
close all force;
set(groot,'defaultFigureVisible',old_visibility);

end


