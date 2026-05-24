% Tests finite-difference and spectral differentiation helpers. Syntax:
%
%                    result=test_finite_difference_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks finite-difference weights, finite-difference matrices,
% Fourier differentiation, Laplacians, FFT differentiation kernels,
% pseudomodulation, and matrix-exponential directional derivatives
% against exact simple cases.
%
% ilya.kuprov@weizmann.ac.il

function result=test_finite_difference_suite()

% Announce the test target
fprintf('TESTING: Finite-difference and spectral differentiation functions\n');

% State the differentiation target of the test
result=new_test_result('kernel/finite_difference_suite',...
                       'Finite-difference and spectral differentiation functions',...
                       'differentiation and pseudomodulation helpers must reproduce exact limiting cases.');

% Three-point centred finite-difference weights at zero are exact and familiar
w=fdweights(0,[-1 0 1],2);
result=test_close(result,'fdweights interpolation',w(1,:),[0 1 0],1e-14,1e-14,...
                  'zeroth-order finite-difference weights interpolate the centre point');
result=test_close(result,'fdweights first derivative',w(2,:),[-1/2 0 1/2],1e-14,1e-14,...
                  'three-point centred first derivative is [-1/2,0,1/2]');
result=test_close(result,'fdweights second derivative',w(3,:),[1 -2 1],1e-14,1e-14,...
                  'three-point centred second derivative is [1,-2,1]');

% Five-point wall finite-difference matrix differentiates quadratics exactly on a unit grid
x=(1:7)'; f=x.^2;
D=fdmat(7,5,1,'wall');
result=test_close(result,'fdmat wall quadratic derivative',D*f,2*x,1e-12,1e-12,...
                  'finite-difference matrix with sufficient stencil differentiates x^2 exactly on a unit grid');
result=test_close(result,'fdvec quadratic derivative',fdvec(f,5,1),2*x,1e-12,1e-12,...
                  'fdvec applies the same finite-difference derivative to a vector');

% Savitzky-Golay differentiation recovers a cubic exactly on a uniform grid
xg=(-5:5)';
f=2*xg.^3-xg.^2+3*xg-1;
df=6*xg.^2-2*xg+3;
sg_der=sgolaydiff(f,1,7,3);
result=test_close(result,'sgolaydiff cubic derivative',sg_der,df,1e-12,1e-12,...
                  'local cubic least-squares differentiation must recover the exact cubic derivative');

% Matrix signal columns must be processed independently
sg_mat=sgolaydiff([f 2*f],1,7,3);
result=test_close(result,'sgolaydiff matrix derivative',sg_mat,[df 2*df],1e-12,1e-12,...
                  'local least-squares differentiation must process matrix columns independently');

% Savitzky-Golay differentiation requires samples down the rows
try
    sgolaydiff(f.',1,7,3);
    error('FAILED: sgolaydiff() accepted a row-vector signal.');
catch err
    result=test_true(result,'sgolaydiff row rejection',...
                     contains(err.message,'sample rows'),...
                     'sgolaydiff() must reject inputs without enough sample rows');
end

% Savitzky-Golay differentiation requires an odd window length
try
    sgolaydiff(f,1,6,3);
    error('FAILED: sgolaydiff() accepted an even window length.');
catch err
    result=test_true(result,'sgolaydiff odd window',...
                     contains(err.message,'npoints must be an odd integer'),...
                     'sgolaydiff() must reject even local least-squares windows');
end

% Savitzky-Golay differentiation suppresses deterministic high-frequency noise
grid=(0:100)'*(2*pi/100);
clean=sin(grid);
noisy=clean+0.03*sin(37*grid)+0.02*cos(29*grid);
raw_der=fdvec(noisy,3,1)/(grid(2)-grid(1));
sg_der=sgolaydiff(noisy,1,15,3)/(grid(2)-grid(1));
raw_err=norm(raw_der-cos(grid),2);
sg_err=norm(sg_der-cos(grid),2);
result=test_true(result,'sgolaydiff noise suppression',sg_err<raw_err,...
                 'local least-squares differentiation should reduce high-frequency noise amplification');

% Periodic finite-difference and Laplacian matrices annihilate constants
Dp=fdmat(8,5,1,'pbc');
result=test_close(result,'fdmat pbc constant derivative',Dp*ones(8,1),zeros(8,1),1e-14,1e-14,...
                  'periodic first derivative of a constant is zero');
Lfd=fdlap([5 4],[1 2],3);
result=test_close(result,'fdlap constant',Lfd*ones(20,1),zeros(20,1),1e-12,1e-12,...
                  'periodic finite-difference Laplacian annihilates constants in multiple dimensions');
Kfd=fdkup([4 4 4],[4 4 4],eye(3),3);
L3=fdlap([4 4 4],[4 4 4],3);
result=test_close(result,'fdkup isotropic identity',Kfd,-L3/3,1e-13,1e-13,...
                  'with an isotropic tensor, fdkup is minus one third of the finite-difference Laplacian');

% Fourier spectral differentiation is exact for a represented sine wave
N=16; [grid,D1]=fourdif(N,1); [~,D2]=fourdif(N,2);
s=sin(grid);
result=test_close(result,'fourdif first derivative',D1*s,cos(grid),1e-12,1e-12,...
                  'Fourier differentiation exactly differentiates sin(x) to cos(x) on the spectral grid');
result=test_close(result,'fourdif second derivative',D2*s,-sin(grid),1e-12,1e-12,...
                  'second Fourier differentiation exactly differentiates sin(x) to -sin(x)');
Lf=fourlap(N,2*pi);
result=test_close(result,'fourlap sine eigenfunction',Lf*s,-s,1e-12,1e-12,...
                  'sin(x) is a Laplacian eigenfunction with eigenvalue -1 on a 2*pi periodic interval');

% FFT differentiation kernel must reproduce the same spectral derivative
kern=fftdiff(1,N,2*pi/N); kern=kern(:);
result=test_close(result,'fftdiff sine derivative',real(ifft(fft(s).*kern)),cos(grid),1e-12,1e-12,...
                  'fftdiff kernel differentiates periodic signals through Fourier multipliers');

% Pseudomodulation zeroth harmonic with zero amplitude is the input spectrum
pm_field=(0:31)'*(2*pi/32);
pm_spec=[sin(2*pm_field) cos(3*pm_field)];
pm_zero=pseudomodulation(pm_field,pm_spec,0,0);
result=test_close(result,'pseudomodulation zero amplitude',pm_zero,pm_spec,1e-14,1e-14,...
                  'zero-amplitude zeroth harmonic must return the input spectrum');

% Pseudomodulation first harmonic tends to half the modulation amplitude times the derivative
pm_amp=1e-4;
pm_first=pseudomodulation(pm_field,sin(2*pm_field),pm_amp,1);
pm_first_ref=pm_amp*cos(2*pm_field);
result=test_close(result,'pseudomodulation first derivative limit',pm_first,pm_first_ref,1e-9,1e-12,...
                  'small-amplitude first harmonic must match the Hyde derivative limit');

% Pseudomodulation second harmonic tends to minus one sixteenth of amplitude squared times the second derivative
pm_second=pseudomodulation(pm_field,sin(2*pm_field),pm_amp,2);
pm_second_ref=(pm_amp^2/4)*sin(2*pm_field);
result=test_close(result,'pseudomodulation second derivative limit',pm_second,pm_second_ref,1e-13,1e-14,...
                  'small-amplitude second harmonic must match the Hyde derivative limit');

% Pseudomodulation requires spectra down the rows
try
    pseudomodulation(pm_field,pm_spec.',pm_amp,1);
    error('FAILED: pseudomodulation() accepted spectra across columns.');
catch err
    result=test_true(result,'pseudomodulation row convention',...
                     contains(err.message,'same number of rows'),...
                     'pseudomodulation() must follow the sgolaydiff() row convention');
end

% Pseudomodulation requires a column field axis
try
    pseudomodulation(pm_field.',pm_spec,pm_amp,1);
    error('FAILED: pseudomodulation() accepted a row-vector field axis.');
catch err
    result=test_true(result,'pseudomodulation field convention',...
                     contains(err.message,'column vector'),...
                     'pseudomodulation() must require a column field axis');
end

% Directional derivative of a commuting matrix exponential has a closed form
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.small_matrix=10;
spin_system.tols.prop_chop=1e-14;
A=diag([1 2]); B=diag([3 4]); Tstep=0.125;
D=dirdiff(spin_system,A,B,Tstep,2);
P=expm(-1i*A*Tstep);
dP=(-1i*Tstep)*B*P;
result=test_close(result,'dirdiff propagator',D{1},P,1e-14,1e-14,...
                  'zeroth directional derivative is the unperturbed propagator');
result=test_close(result,'dirdiff commuting derivative',D{2},dP,1e-14,1e-14,...
                  'for commuting A and B, the first derivative is -i*T*B*exp(-i*A*T)');

end
