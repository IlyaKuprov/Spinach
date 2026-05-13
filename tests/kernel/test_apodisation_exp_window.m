% Tests exponential FID apodisation. Syntax:
%
%                    result=test_apodisation_exp_window()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the explicit exponential window and the NMR Fourier
% convention that halves the first point of each active FID dimension.
%
% ilya.kuprov@weizmann.ac.il

function result=test_apodisation_exp_window()

% State the processing target of the test
result=new_test_result('kernel/apodisation_exp_window',...
                       'Exponential FID apodisation',...
                       'exponential apodisation must multiply by exp(-k*x) and halve the first point.');

% Build a minimal reporting object and a constant FID
spin_system.sys.output='hush';
fid=ones(4,1);

% Apply an exponential window
fid_obs=apodisation(spin_system,fid,{{'exp',1}});
fid_ref=exp(-linspace(0,1,4)).';
fid_ref(1)=fid_ref(1)/2;

% Check the explicit window
result=test_close(result,'exp window and first-point half',fid_obs,fid_ref,1e-15,1e-15,...
                  'the first point is halved, then multiplied by exp(-x)');

end

