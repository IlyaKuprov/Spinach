% Tests scalar unit-conversion functions. Syntax:
%
%                    result=test_unit_conversion_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks physically defined unit conversions used across magnetic
% resonance: Hartree to J/mol, cm^-1 to Hz, field/frequency conversions,
% and Lorentzian linewidth to R2.
%
% ilya.kuprov@weizmann.ac.il

function result=test_unit_conversion_suite()

% State the conversion target of the test
result=new_test_result('kernel/unit_conversion_suite',...
                       'Unit-conversion functions',...
                       'unit conversion helpers must implement their defining physical constants.');

% Check Hartree energy conversion
hartree=[0 1 2.5];
result=test_close(result,'hartree2joule',hartree2joule(hartree),2625499.62*hartree,1e-10,1e-15,...
                  'one Hartree is 2625499.62 J/mol in Spinach convention');

% Check inverse-centimetre and frequency conversion
icm=[0 1 12.5];
hz=100*299792458*icm;
result=test_close(result,'icm2hz',icm2hz(icm),hz,1e-6,1e-15,...
                  'one inverse centimetre is c*100 Hz');
result=test_close(result,'hz2icm',hz2icm(hz),icm,1e-12,1e-12,...
                  'hz2icm is the algebraic inverse of icm2hz');

% Check electron-field hyperfine conversions
hfc_gauss=[0 10 25];
g=2.0023193043622;
muB=9.274009994e-24;
hbar=1.054571628e-34;
conv=1e-10*g*muB/(hbar*2*pi);
hfc_mhz=conv*hfc_gauss;
result=test_close(result,'gauss2mhz',gauss2mhz(hfc_gauss,g),hfc_mhz,1e-12,1e-12,...
                  'Gauss hyperfine units are converted through the electron Zeeman frequency');
result=test_close(result,'mhz2gauss',mhz2gauss(hfc_mhz,g),hfc_gauss,1e-12,1e-12,...
                  'mhz2gauss is the algebraic inverse of gauss2mhz');

% Check milliTesla to Hz conversion
hfc_mt=[0 1 3.5];
hfc_hz=1e-3*g*muB*hfc_mt/(hbar*2*pi);
result=test_close(result,'mt2hz',mt2hz(hfc_mt,g),hfc_hz,1e-6,1e-12,...
                  'milliTesla hyperfine units are converted to linear frequency through g*muB/hbar');

% Check Lorentzian linewidth conversion
fwhm=[1 2.5 10];
result=test_close(result,'fwhm2rlx',fwhm2rlx(fwhm),pi*fwhm,1e-15,1e-15,...
                  'Lorentzian full width at half maximum corresponds to R2=pi*FWHM');

end

