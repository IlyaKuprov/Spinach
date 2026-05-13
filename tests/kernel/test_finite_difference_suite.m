% Tests finite-difference and spectral differentiation helpers. Syntax:
%
%                    result=test_finite_difference_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks finite-difference weights, finite-difference matrices,
% Fourier differentiation, Laplacians, FFT differentiation kernels, and
% matrix-exponential directional derivatives against exact simple cases.
%
% ilya.kuprov@weizmann.ac.il

function result=test_finite_difference_suite()

% Announce the test target
fprintf('TESTING: Finite-difference and spectral differentiation functions\n');

% State the differentiation target of the test
result=new_test_result('kernel/finite_difference_suite',...
                       'Finite-difference and spectral differentiation functions',...
                       'differentiation helpers must reproduce exact low-order polynomial and Fourier derivatives.');

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
