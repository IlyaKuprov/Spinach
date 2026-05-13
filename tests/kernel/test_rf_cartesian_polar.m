% Tests RF Cartesian and polar waveform conversion. Syntax:
%
%                    result=test_rf_cartesian_polar()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that RF amplitude/phase coordinates round-trip to X/Y
% controls and that gradients transform by the chain rule.
%
% ilya.kuprov@weizmann.ac.il

function result=test_rf_cartesian_polar()

% Announce the test target
fprintf('TESTING: RF Cartesian/polar conversion\n');

% State the pulse-control target of the test
result=new_test_result('kernel/rf_cartesian_polar',...
                       'RF Cartesian/polar conversion',...
                       'Cartesian and polar RF controls must describe the same complex waveform.');

% Define a waveform away from the zero-amplitude singularity
r=[1.0 2.0 3.0];
p=[0.0 pi/3 -pi/2];
Dr=[4.0 5.0 6.0];
Dp=[0.5 -0.25 0.75];

% Convert to Cartesian and back
[x,y,Dx,Dy]=polar2cartesian(r,p,Dr,Dp);
[r_back,p_back,Dr_back,Dp_back]=cartesian2polar(x,y,Dx,Dy);

% Check coordinate round-trip
result=test_close(result,'RF x coordinate',x,r.*cos(p),1e-15,1e-15,...
                  'x control is amplitude times cosine phase');
result=test_close(result,'RF y coordinate',y,r.*sin(p),1e-15,1e-15,...
                  'y control is amplitude times sine phase');
result=test_close(result,'amplitude round-trip',r_back,r,1e-15,1e-15,...
                  'non-zero RF amplitude is invariant under coordinate round-trip');
result=test_close(result,'phase round-trip',p_back,p,1e-15,1e-15,...
                  'phase is recovered by atan2 for the same quadrant');

% Check gradient round-trip
result=test_close(result,'amplitude gradient round-trip',Dr_back,Dr,1e-12,1e-12,...
                  'gradient components transform by the chain rule');
result=test_close(result,'phase gradient round-trip',Dp_back,Dp,1e-12,1e-12,...
                  'phase derivative is the angular component of the same gradient');

end

