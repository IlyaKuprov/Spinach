% FSLG-HETCOR of alpha-glycine powder under MAS.
%
% Calculation time: hours on NVidia Tesla A100,
%                   much longer on CPU
%
% guinevere.mathies@uni-konstanz.de

function fslghetcor_mas_gly()

% Spin system properties (PCM DFT calculation)
[sys,inter]=g2spinach(gparse('..\..\examples\standard_systems\glycine.log'),...
                      {{'H','1H'},{'C','13C'}},[31.8 182.1],[]);

% 400 MHz spectrometer
sys.magnet=9.4;

% Isotropic alpha-glycine chemical shifts
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4); % CO
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,43.6);  % CA
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,3,2.6);   % H_CA
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,4,3.8);   % H_CA
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,5,8.0);   % H_N
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,6,8.0);   % H_N
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,7,8.0);   % H_N

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=4;

% Ignore interactions below 200 Hz
sys.tols.inter_cutoff=2*pi*200;

% Use GPU arithmetic
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Start with Lz of 1H, detect in quadrature on 13C
parameters.rho0=state(spin_system,'Lz','1H');
parameters.coil=state(spin_system,'L+','13C');

% Experiment setup
parameters.spins={'1H','13C'};
parameters.rate=10000;
parameters.axis=[1 1 1];
parameters.max_rank=7;
parameters.grid='rep_2ang_100pts_sph';
parameters.offset=[2e3 1e4];      % Hz
parameters.hi_pwr=83e3;           % Hz
parameters.cp_pwr=[60e3 50e3];    % Hz
parameters.cp_dur=1e-4;           % seconds
parameters.sweep=[NaN 1/(33e-6)]; % f1 f2
parameters.npoints=[128 512];     % f1 f2
parameters.zerofill=[512 2048];   % f1 f2
parameters.nblocks=4;             % Number of FSLG blocks
parameters.verbose=1;
parameters.invert_axis=1;

% Simulation
fid=singlerot(spin_system,@fslghetcor,parameters,'nmr');

% Apodisation
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1);
f1_sin=fftshift(fft(fid.sin,parameters.zerofill(2),1),1);

% States quadrature and F1 Fourier transform
f1_states=real(f1_cos)+1i*real(f1_sin);
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Custom sweep and offset in the F1 direction: 2 pulses per block 
% in frequency-switched Lee-Goldburg; 1/sqrt(3) is the scaling of 
% the chemical shift in the effective Hamiltonian, it may vary in
% the actual experiment; sqrt(2/3) accounts for the fact that the
% effective field is at the magic angle.
parameters.sweep(1)=sqrt(3)*parameters.hi_pwr/(2*parameters.nblocks*sqrt(2/3));
parameters.offset(1)=0;

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 0.95 0.05 0.95],2,256,6,'both');

end

