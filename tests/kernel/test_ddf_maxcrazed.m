% Check the two-dimensional MAXCRAZED example against the qualitative
% Fig. 1(c) / Table I discrimination. Syntax:
%
%                    result=test_ddf_maxcrazed()
%
% Outputs:
%
%   result - regression test result for peak assignment and contrast

function result=test_ddf_maxcrazed()

% Run the production example with the caller's original path restored on exit
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_path=path;
restore_path=onCleanup(@()path(old_path));
addpath(fullfile(root_dir,'examples','imaging','distant_dipolar'));
data=ddf_maxcrazed();

% Verify the four double-quantum and single-quantum frequency assignments
result=new_test_result('kernel/ddf_maxcrazed','Two-axis MAXCRAZED',...
                       'The four spectral peaks must distinguish XY and XZ encodings.');
[~,strongest]=maxk(reshape(data.spectrum(:,:,1),[],1),4);
[f1,f2]=ind2sub([numel(data.f1),numel(data.f2)],strongest);
found=sortrows([reshape(data.f1(f1),[],1) reshape(data.f2(f2),[],1)]);
expected=[1260  630;...
          2130  630;...
          2130 1500;...
          3000 1500];
result=test_close(result,'four X-CRAZED spectral peaks',found,expected,1e-10,0,...
                  'the four strongest peaks resolve at the Table I frequency pairs');

% Check the anisotropy without asserting a fit to the measured heights
predicted=data.peaks.predicted;
result=test_true(result,'both self peaks are suppressed',...
                 all(predicted([1 4],3)<predicted([1 4],2)),...
                 'the strong Z gradient suppresses the self torque relative to Y');
result=test_true(result,'both mixed peaks are retained',...
                 all(predicted(2:3,3)>=predicted(2:3,2)),...
                 'the strong Z gradient does not attenuate the mixed peaks relative to Y');
result=test_true(result,'hetero/homo ratio more than doubles',...
                 data.predicted_ratio(3)>2*data.predicted_ratio(2),...
                 'the change in gradient orientation is not simple diffusion attenuation');

end


