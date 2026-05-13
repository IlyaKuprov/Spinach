% Tests dynamic pulse utility paths. Syntax:
%
%               result=test_dynamic_pulse_utilities_deep()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks gradient-pulse dynamics, heterodyne filtering, RLC
% response transforms, Bruker pulse-file writing, finite-RF R-sequence
% compilation, waveform basis variants, and pulse-shape variants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_pulse_utilities_deep()

% State the pulse-utility target of the test
result=new_test_result('kernel/dynamic_pulse_utilities_deep',...
                       'Dynamic pulse utility paths',...
                       'pulse utilities must produce finite dynamic outputs and match direct small-matrix references.');

% Build a one-proton Liouville-space spin system with a carrier frequency
spin_system=local_liouv_system(1.0);
L=0*operator(spin_system,'Lz',1);
rho=state(spin_system,'Lx',1);

% Compare grad_pulse with a direct small-matrix exponential reference
rho_grad=grad_pulse(spin_system,L,rho,4.0,0.12,1.5e-4,0.8);
rho_ref=grad_pulse_ref(spin_system,L,rho,4.0,0.12,1.5e-4,0.8);
result=test_close(result,'grad_pulse direct exponential reference',rho_grad,rho_ref,1e-10,1e-10,...
                  'Edwards single-gradient propagation is the corresponding auxiliary-matrix exponential');
result=test_true(result,'grad_pulse non-zero response',norm(rho_grad-rho,2)>1e-8,...
                 'a finite gradient changes transverse magnetisation in the one-spin carrier frame');

% Compare grad_sandw with a direct small-matrix exponential reference
pulse_prop=propagator(spin_system,2*pi*50*operator(spin_system,'Ly',1),2.0e-4);
rho_sandw=grad_sandw(spin_system,L,rho,pulse_prop,[3.0 -2.0],0.10,[1.0e-4 1.4e-4],[0.7 0.9]);
rho_ref=grad_sandw_ref(spin_system,L,rho,pulse_prop,[3.0 -2.0],0.10,[1.0e-4 1.4e-4],[0.7 0.9]);
result=test_close(result,'grad_sandw direct exponential reference',rho_sandw,rho_ref,1e-10,1e-10,...
                  'Edwards gradient-sandwich propagation is the corresponding block auxiliary exponential');
result=test_true(result,'grad_sandw non-zero response',norm(rho_sandw-pulse_prop*rho,2)>1e-8,...
                 'finite sandwich gradients change the propagated transverse magnetisation');

% Heterodyne a clean wall-clock carrier into the rotating frame
dt=1e-4;
freq=1000;
time_grid=dt*((1:4096)'-1);
signal=cos(2*pi*freq*time_grid);
[X_het,Y_het]=heterodyne(dt,signal,freq);
steady=200:(numel(signal)-200);
result=test_close(result,'heterodyne in-phase mean',mean(X_het(steady)),1,5e-2,5e-2,...
                  'mixing a unit cosine with its carrier leaves unit in-phase DC signal after filtering');
result=test_close(result,'heterodyne quadrature mean',mean(Y_het(steady)),0,5e-2,5e-2,...
                  'a pure cosine has zero quadrature DC component in the rotating frame');
result=test_true(result,'heterodyne finite output',all(isfinite([X_het;Y_het])),...
                 'heterodyne output samples must be finite real numbers');

% Exercise piecewise-constant and piecewise-linear RLC response models
omega=2*pi*800;
Q=8;
dt_user=1e-3;
[X_pwc,Y_pwc,dt_pwc]=restrans([1;2;1;0],[0;0.5;0;-0.5],dt_user,omega,Q,'pwc',2);
[X_pwl,Y_pwl,dt_pwl]=restrans([1;2;1;0;-1],[0;0.5;0;-0.5;0],dt_user,omega,Q,'pwl',2);
[X_tsc,Y_tsc,dt_tsc]=restrans([1;2;1;0;-1],[0;0.5;0;-0.5;0],dt_user,omega,Q,'pwl_tsc',2);
result=test_true(result,'restrans pwc finite output',iscolumn(X_pwc)&&iscolumn(Y_pwc)&&...
                 (dt_pwc>0)&&all(isfinite([X_pwc;Y_pwc])),...
                 'piecewise-constant RLC response returns finite column waveforms and a positive time step');
result=test_true(result,'restrans pwl finite output',iscolumn(X_pwl)&&iscolumn(Y_pwl)&&...
                 (dt_pwl>0)&&all(isfinite([X_pwl;Y_pwl])),...
                 'piecewise-linear RLC response returns finite column waveforms and a positive time step');
result=test_close(result,'restrans pwl_tsc samples',X_tsc+1i*Y_tsc,X_pwl+1i*Y_pwl,1e-12,1e-12,...
                  'time-shift compensation changes the diagnostic time grid but not returned waveform samples');
result=test_close(result,'restrans pwl_tsc timestep',dt_tsc,dt_pwl,1e-15,1e-15,...
                  'time-shift compensation leaves the returned waveform time step unchanged');

% Write a temporary Bruker pulse file and inspect stable JCAMP fields
file_name=[tempname(fileparts(mfilename('fullpath'))) '.txt'];
cleanup=onCleanup(@()local_delete(file_name));
bruker_write([1;0;-1],[0;1;0],2.5e-6,file_name);
lines=readlines(file_name);
result=test_true(result,'bruker_write file exists',exist(file_name,'file')==2,...
                 'bruker_write must create the requested temporary output file');
result=test_true(result,'bruker_write npoints',any(contains(lines,"##NPOINTS=3")),...
                 'the Bruker file declares the number of waveform points');
result=test_true(result,'bruker_write shape length',any(contains(lines,"##$SHAPE_LENGTH=7.5")),...
                 'the Bruker file declares the total pulse length in microseconds');
result=test_true(result,'bruker_write terminator',any(strcmp(strtrim(lines),"##END")),...
                 'the Bruker file is terminated by the JCAMP end marker');

% Check finite-RF R-sequence compiler branches against direct exponentials
spin_system_h=local_hilb_system();
S=pauli(2);
Sx=S.x;
Sy=S.y;
pulse_amp=7.5;
pulse_dur=0.03;
[P,T]=rseq_compiler(spin_system_h,zeros(2),Sx,Sy,[0;pi/2;0],pulse_amp,pulse_dur,'180_pulse');
result=test_close(result,'rseq_compiler finite index map',T,[1;2;1],1e-15,1e-15,...
                  'finite-RF compilation still reuses propagators for repeated phases');
result=test_close(result,'rseq_compiler finite x pulse',P{1},expm(-1i*pulse_amp*Sx*pulse_dur),1e-14,1e-14,...
                  'a zero-phase finite-RF element exponentiates the X pulse Hamiltonian');
result=test_close(result,'rseq_compiler finite y pulse',P{2},expm(-1i*pulse_amp*Sy*pulse_dur),1e-14,1e-14,...
                  'a pi/2 finite-RF element exponentiates the Y pulse Hamiltonian');
[P,T]=rseq_compiler(spin_system_h,zeros(2),Sx,Sy,[0;pi/2],pulse_amp,[0.01 0.02],'90270_pulse');
result=test_close(result,'rseq_compiler composite index map',T,[1;2],1e-15,1e-15,...
                  'composite finite-RF elements compile the two unique phases separately');
result=test_close(result,'rseq_compiler composite x pulse',P{1},expm(-1i*pulse_amp*Sx*0.01),1e-14,1e-14,...
                  'the first composite finite-RF pulse uses the first duration');
result=test_close(result,'rseq_compiler composite y pulse',P{2},expm(-1i*pulse_amp*Sy*0.02),1e-14,1e-14,...
                  'the second composite finite-RF pulse uses the second duration');

% Check waveform basis variants on a compact odd grid
basis_types={'sine_waves','cosine_waves','legendre'};
for n=1:numel(basis_types)

    % Build and check each supported basis family
    basis_waves=wave_basis(basis_types{n},3,17);
    result=test_close(result,[basis_types{n} ' deep Gram matrix'],basis_waves'*basis_waves,eye(3),1e-12,1e-12,...
                      'supported waveform basis variants return orthonormal columns');
    result=test_true(result,[basis_types{n} ' deep dimensions'],isequal(size(basis_waves),[17 3]),...
                     'basis matrices must have one row per time point and one column per function');
end

% Check all supported pulse-shape variants against explicit formulae
result=test_close(result,'pulse_shape rectangular deep',pulse_shape('rectangular',5),ones(1,5),1e-15,1e-15,...
                  'rectangular pulse-shape samples are all unity');
result=test_close(result,'pulse_shape sinc3 deep',pulse_shape('sinc3',5),pi*sinc(linspace(-3,3,5)),1e-15,1e-15,...
                  'sinc3 pulse-shape samples follow pi*sinc over [-3,3]');
result=test_close(result,'pulse_shape sinc5 deep',pulse_shape('sinc5',5),pi*sinc(linspace(-5,5,5)),1e-15,1e-15,...
                  'sinc5 pulse-shape samples follow pi*sinc over [-5,5]');
time_grid=linspace(-2,2,5);
gauss_ref=exp(-(time_grid.^2)/2)/(sqrt(2*pi)*sqrt(2));
result=test_close(result,'pulse_shape gaussian deep',pulse_shape('gaussian',5),gauss_ref,1e-15,1e-15,...
                  'Gaussian pulse-shape samples follow normpdf(t)/sqrt(2) over [-2,2]');

end

% Local quiet one-spin Liouville-space test system
function spin_system=local_liouv_system(magnet)

% Specify the spin system
sys.magnet=magnet;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};

% Specify a full Liouville basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Build the quiet regression-test system
spin_system=test_spin_system(sys,inter,bas);

end

% Local quiet one-spin Hilbert-space test system
function spin_system=local_hilb_system()

% Specify the spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};

% Specify a full Hilbert basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Build the quiet regression-test system
spin_system=test_spin_system(sys,inter,bas);

end

% Direct small-matrix reference for a single gradient pulse
function rho=grad_pulse_ref(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)

% Build the effective gradient operator
G=1e-4*s_fac*g_amp*s_len*g_dur*carrier(spin_system,'all')/spin_system.inter.magnet;

% Propagate through the single-gradient auxiliary system
rho=expm(-1i*full(L)*g_dur)*rho;
rho=expm(1i*full(G)/2)*rho;
aux_mat=[0*G,1i*speye(size(G));0*G,G];
aux_rho=expm(-1i*full(aux_mat))*[0*rho;rho];
rho=aux_rho(1:end/2,:);

end

% Direct small-matrix reference for a gradient sandwich
function rho=grad_sandw_ref(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)

% Build the effective gradient operators
R=carrier(spin_system,'all')/spin_system.inter.magnet;
G1=1e-4*s_facs(1)*g_amps(1)*s_len*g_durs(1)*R;
G2=1e-4*s_facs(2)*g_amps(2)*s_len*g_durs(2)*R;

% Propagate through the gradient-sandwich auxiliary system
rho=expm(-1i*full(L)*g_durs(1))*rho;
rho=expm(1i*full(G1)/2)*rho;
aux_mat=[-G2,1i*P;0*G1,G1];
aux_rho=expm(-1i*full(aux_mat))*[0*rho;rho];
rho=aux_rho(1:end/2,:);
rho=expm(-1i*full(G2)/2)*rho;
rho=expm(-1i*full(L)*g_durs(2))*rho;

end

% Best-effort temporary file deletion
function local_delete(file_name)

% Delete the file if it exists
if exist(file_name,'file')==2
    delete(file_name);
end

end


