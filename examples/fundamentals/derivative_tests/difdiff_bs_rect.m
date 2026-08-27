% Directional derivative test for Cartesian GRAPE with Bloch-Siegert
% corrections. 
%
% aditya.dev@weizmann.ac.il

function difdiff_bs_rect()

% Set the magnetic field
sys.magnet=10.2;

% Set isotopes
sys.isotopes={'1H','1H','13C','13C'};

% Set interactions
inter.zeeman.scalar={1.5,2.0,30.0,40.0};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=7.0;
inter.coupling.scalar{1,3}=150;
inter.coupling.scalar{2,4}=150;
inter.coupling.scalar{3,4}=50;

% Set basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Build and normalise initial state
rho_init=singlet(spin_system,1,2);
rho_init=rho_init/norm(full(rho_init),2);

% Build and normalise target state
rho_targ=state(spin_system,{'Lz'},{4});
rho_targ=rho_targ/norm(full(rho_targ),2);

% Get control operators
LxH=operator(spin_system,'Lx','1H');
LyH=operator(spin_system,'Ly','1H');
LxC=operator(spin_system,'Lx','13C');
LyC=operator(spin_system,'Ly','13C');

% Get offset and shift operators
LzH=operator(spin_system,'Lz','1H');
LzC=operator(spin_system,'Lz','13C');

% Build drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control structure
control.drifts={{H}};
control.operators={LxH,LyH,LxC,LyC};
control.rho_init={rho_init};
control.rho_targ={rho_targ};
control.pwr_levels=2*pi*500;
control.max_iter=1000;
control.plotting={};
control.integrator='rectangle';
control.pulse_dt=1.5e-4*ones(1,100);

% Define transmitter offsets for both isotopes
control.offsets={1050,5285};
control.off_ops={LzH,LzC};

% Enable Bloch-Siegert corrections
control.bsiegert=true();
control.isotopes={'1H','13C'};
control.channels=[1,1,2,2];

% Run Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Build random guess and finite difference increment
guess=randn(numel(control.operators),numel(control.pulse_dt))/10;
h=1e-5;
fail_count=0;

% Request analytical gradient
[~,~,grad_anl]=grape_xy(guess,spin_system);
grad_anl=squeeze(grad_anl(:,:,1));

% Left waveform edge
wave_forw=guess;
wave_forw(1)=wave_forw(1)+h;
wave_back=guess;
wave_back(1)=wave_back(1)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
rel_err=abs(grad_anl(1)-grad_num)/max(abs(grad_num),eps('double'));
if rel_err<5e-6
    disp('left edge test passed');
else
    disp(['left edge test failed, grad_anl=' num2str(grad_anl(1),'%.12e') ...
        ', grad_num=' num2str(grad_num,'%.12e') ', rel_err=' ...
        num2str(rel_err,'%.12e')]);
    fail_count=fail_count+1;
end

% Right waveform edge
wave_forw=guess;
wave_forw(end)=wave_forw(end)+h;
wave_back=guess;
wave_back(end)=wave_back(end)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
rel_err=abs(grad_anl(end)-grad_num)/max(abs(grad_num),eps('double'));
if rel_err<5e-6
    disp('right edge test passed');
else
    disp(['right edge test failed, grad_anl=' num2str(grad_anl(end),'%.12e') ...
        ', grad_num=' num2str(grad_num,'%.12e') ', rel_err=' ...
        num2str(rel_err,'%.12e')]);
    fail_count=fail_count+1;
end

% Waveform midpoint
n_rows=size(guess,1);
n_cols=size(guess,2);
mid_row=ceil(n_rows/2);
mid_col=ceil(n_cols/2);
mid_idx=sub2ind([n_rows n_cols],mid_row,mid_col);
wave_forw=guess;
wave_forw(mid_idx)=wave_forw(mid_idx)+h;
wave_back=guess;
wave_back(mid_idx)=wave_back(mid_idx)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
rel_err=abs(grad_anl(mid_idx)-grad_num)/max(abs(grad_num),eps('double'));
if rel_err<5e-6
    disp('midpoint test passed');
else
    disp(['midpoint test failed, grad_anl=' num2str(grad_anl(mid_idx),'%.12e') ...
        ', grad_num=' num2str(grad_num,'%.12e') ', rel_err=' ...
        num2str(rel_err,'%.12e')]);
    fail_count=fail_count+1;
end

% Print final summary
if fail_count==0
    disp('all derivative tests passed');
else
    error([int2str(fail_count) ' derivative test(s) failed']);
end

end


