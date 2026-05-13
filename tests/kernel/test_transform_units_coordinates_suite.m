% Tests unit and coordinate transform helpers. Syntax:
%
%                    result=test_transform_units_coordinates_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks scalar physical constants, inverse unit conversions,
% crystallographic coordinate conversion, and ISO spherical coordinates.
%
% ilya.kuprov@weizmann.ac.il

function result=test_transform_units_coordinates_suite()

% State the conversion target of the test
result=new_test_result('kernel/transform_units_coordinates_suite',...
                       'Unit and coordinate transforms',...
                       'unit and coordinate transforms must implement their defining constants and inverse maps.');

% Check Hartree energy conversion
hartree=[0 1 2.5];
result=test_close(result,'hartree2joule constant',hartree2joule(hartree),2625499.62*hartree,1e-10,1e-15,...
                  'one Hartree is 2625499.62 J/mol in the Spinach convention');

% Check inverse-centimetre and Hz conversions
icm=[0 1 12.5];
hz=100*299792458*icm;
result=test_close(result,'icm2hz constant',icm2hz(icm),hz,1e-6,1e-15,...
                  'one inverse centimetre is c*100 Hz');
result=test_close(result,'hz2icm inverse',hz2icm(hz),icm,1e-12,1e-12,...
                  'hz2icm must be the algebraic inverse of icm2hz');

% Check Angstrom^3 and cgs-ppm susceptibility conversion
ang=[-2 0 3.5];
cgsppm=6.02214129e23*ang/(4*pi*1e18);
result=test_close(result,'ang2cgsppm constant',ang2cgsppm(ang),cgsppm,1e-15,1e-15,...
                  'cubic Angstrom susceptibility is converted to cm^3/mol through Avogadro number and 4*pi');
result=test_close(result,'cgsppm2ang inverse',cgsppm2ang(cgsppm),ang,1e-14,1e-14,...
                  'cgsppm2ang must invert ang2cgsppm element by element');

% Check chemical shift and frequency conversion including isotope sign
ppm=[-1 0 3.2];
B0=14.1;
hz=ppm2hz(ppm,B0,'1H');
result=test_close(result,'ppm2hz hz2ppm inverse',hz2ppm(hz,B0,'1H'),ppm,1e-12,1e-12,...
                  'chemical shift offsets must round-trip at fixed field and isotope');

% Check electron field-frequency conversions
hfc_gauss=[0 10 25];
g=2.0023193043622;
muB=9.274009994e-24;
hbar=1.054571628e-34;
conv=1e-10*g*muB/(hbar*2*pi);
hfc_mhz=conv*hfc_gauss;
result=test_close(result,'gauss2mhz constant',gauss2mhz(hfc_gauss,g),hfc_mhz,1e-12,1e-12,...
                  'Gauss hyperfine units are converted through the electron Zeeman frequency in MHz');
result=test_close(result,'mhz2gauss inverse',mhz2gauss(hfc_mhz,g),hfc_gauss,1e-12,1e-12,...
                  'mhz2gauss must invert gauss2mhz for the same g factor');

% Check milliTesla to Hz and g-value to frequency conversion
hfc_mt=[0 1 3.5];
hfc_hz=1e-3*g*muB*hfc_mt/(hbar*2*pi);
result=test_close(result,'mt2hz constant',mt2hz(hfc_mt,g),hfc_hz,1e-6,1e-12,...
                  'milliTesla hyperfine units convert to linear frequency through g*muB/hbar');
gvals=[2.0023193043622 2.1];
B=0.34;
freq_ref=B*spin('E')*gvals/(2*pi*2.0023193043622);
result=test_close(result,'g2freq proportionality',g2freq(gvals,B),freq_ref,1e-6,1e-12,...
                  'g-value frequencies scale linearly with magnetic field and with g/ge');

% Check Lorentzian full-width to transverse relaxation rate conversion
fwhm=[1 2.5 10];
result=test_close(result,'fwhm2rlx constant',fwhm2rlx(fwhm),pi*fwhm,1e-15,1e-15,...
                  'a Lorentzian line with FWHM in Hz has R2=pi*FWHM');

% Check orthorhombic fractional-to-Cartesian crystal coordinates
ABC=[0 0 0;1/2 1/3 1/4;1 1 1];
[XYZ,va,vb,vc]=frac2cart(2,3,4,90,90,90,ABC);
XYZ_ref=[0 0 0;1 1 1;2 3 4];
result=test_close(result,'frac2cart orthorhombic coordinates',XYZ,XYZ_ref,1e-14,1e-14,...
                  'for a 90-degree unit cell fractional coordinates scale by a, b, and c');
result=test_close(result,'frac2cart primitive a',va,[2;0;0],1e-14,1e-14,...
                  'the first primitive vector of an orthorhombic cell is along x');
result=test_close(result,'frac2cart primitive b',vb,[0;3;0],1e-14,1e-14,...
                  'the second primitive vector of an orthorhombic cell is along y');
result=test_close(result,'frac2cart primitive c',vc,[0;0;4],1e-14,1e-14,...
                  'the third primitive vector of an orthorhombic cell is along z');

% Check ISO spherical coordinate convention on Cartesian axes
x=[1;0;0]; y=[0;1;0]; z=[0;0;1];
[r,theta,phi]=xyz2sph(x,y,z);
result=test_close(result,'xyz2sph radii',r,ones(3,1),1e-15,1e-15,...
                  'unit Cartesian basis vectors all have radius one');
result=test_close(result,'xyz2sph inclinations',theta,[pi/2;pi/2;0],1e-15,1e-15,...
                  'the ISO inclination is measured down from the positive z axis');
result=test_close(result,'xyz2sph azimuths',phi,[0;pi/2;0],1e-15,1e-15,...
                  'the ISO azimuth is atan2(y,x) in the xy plane');

end


