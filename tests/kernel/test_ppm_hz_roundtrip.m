% Tests chemical shift and frequency conversion. Syntax:
%
%                    result=test_ppm_hz_roundtrip()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the physical definition nu=delta*1e-6*gamma*B0/(2*pi)
% and verifies that the inverse conversion preserves the sign of gamma.
%
% ilya.kuprov@weizmann.ac.il

function result=test_ppm_hz_roundtrip()

% Announce the test target
fprintf('TESTING: Chemical-shift frequency conversion\n');

% State the physical target of the test
result=new_test_result('kernel/ppm_hz_roundtrip',...
                       'Chemical-shift frequency conversion',...
                       'ppm2hz and hz2ppm must implement the Larmor-frequency definition.');

% Define a field and shifts for positive-gamma and negative-gamma nuclei
B0=14.1;
ppm=[-2.5 0 3.0 12.0];

% Check explicit frequency formula for proton shifts
hz_ref=1e-6*ppm*(B0*spin('1H')/(2*pi));
hz_obs=ppm2hz(ppm,B0,'1H');
result=test_close(result,'1H explicit formula',hz_obs,hz_ref,1e-9,1e-14,...
                  'chemical shift in ppm is fractional Larmor offset');
result=test_close(result,'1H inverse conversion',hz2ppm(hz_obs,B0,'1H'),ppm,1e-12,1e-12,...
                  'frequency-to-ppm conversion is the algebraic inverse');

% Check that negative magnetogyric ratios retain sign
hz_ref=1e-6*ppm*(B0*spin('15N')/(2*pi));
hz_obs=ppm2hz(ppm,B0,'15N');
result=test_close(result,'15N sign preservation',hz_obs,hz_ref,1e-9,1e-14,...
                  'negative-gamma nuclei must produce negative offsets for positive ppm');
result=test_close(result,'15N inverse conversion',hz2ppm(hz_obs,B0,'15N'),ppm,1e-12,1e-12,...
                  'inverse conversion must preserve the sign convention');

end

