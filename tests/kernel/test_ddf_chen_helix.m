% Signed phase-sensitive tests of the Chen transverse-helix limit.
%
%                     result=test_ddf_chen_helix()
%
% Output:
%
%   result - checks of axial shifts, magic-angle nulls, and receiver sign

function result=test_ddf_chen_helix()

% Find the example and restore the caller path on return
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_path=path;
restore_path=onCleanup(@()path(old_path));
addpath(fullfile(root_dir,'examples','imaging','distant_dipolar'));

% Run the production example and compare with the single-helix solution
result=new_test_result('kernel/ddf_chen_helix','Chen transverse helix',...
                       'Spherical water with weak and inverted transverse helices.');
data=ddf_chen_helix();
target=-cosd(data.tips')/(3*0.130);
result=test_close(result,'signed axial slopes',data.omega(:,1),target,...
                  1e-4,0,'the transverse field torques the uniform longitudinal moment');

% Check the full complex echo against rotation at the predicted rate
echo=exp(1i*data.time*target.');
result=test_close(result,'axial phase and amplitude',max(abs(data.echo(:,:,1)-echo),[],'all'),...
                  0,1e-4,0,'one spatial mode and a constant longitudinal moment');
result=test_close(result,'magic-angle null',max(abs(data.echo(:,:,2)-1),[],'all'),...
                  0,1e-8,0,'both the uniform and helical field multipliers vanish');
result=test_close(result,'offset sign',data.offset_omega,-2*pi,...
                  1e-7,0,'positive Spinach offset rotates L+ detection negatively');

end


