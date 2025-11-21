% A soft pulse simulation for a gadolinium ion. The soft pulse
% is simulated using Fokker-Planck formalism. Zero-field split-
% ting distribution is sampled using the statistical parameters
% reported in Figure 5 of Raitsimring et al, App. Mag. Res. 28,
% 281-295 (2005). Powder average simulation with a third-order 
% numerical rotating frame transformation.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function spa_gd_dota_powder()

% Preallocate the spectrum
spectrum=zeros(2048,1,'like',1i);

% Get the sampling
[D,E,W]=zfs_sampling(30,5,1e-2); drawnow;

% Get the figure going
kfigure();

% Loop over ZFS distribution
for n=1:numel(W)
    
    % Spin system parameters
    sys.magnet=3.5;
    sys.isotopes={'E8'};
    inter.zeeman.scalar={2.002319};
    inter.coupling.matrix{1,1}=0.56e9*zfs2mat(D(n),E(n),0,0,0);
    
    % Basis set
    bas.formalism='sphten-liouv';
    bas.approximation='none';
    bas.projections=-3:3;
    
    % Disable trajectory-level SSR algorithms
    sys.disable={'trajlevel'};
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Sequence parameters
    parameters.spins={'E8'};
    parameters.rho0=state(spin_system,'Lz','E8');
    parameters.coil=state(spin_system,'L+','E8');
    parameters.decouple={};
    parameters.offset=0;
    parameters.sweep=0.8e10;
    parameters.npoints=512;
    parameters.zerofill=2048;
    parameters.axis_units='GHz';
    parameters.grid='rep_2ang_400pts_sph';
    parameters.derivative=0;
    parameters.invert_axis=0;
    parameters.rframes={{'E8',3}};
    parameters.verbose=0;
    
    % Soft pulse parameters
    parameters.pulse_rnk=2;
    parameters.pulse_phi=-pi/2;
    parameters.pulse_frq=-0.5e9;
    parameters.pulse_dur=50.0e-9;
    parameters.pulse_pwr=2*pi*0.02e+9;
    parameters.method='expm';
    
    % Simulation 
    fid=powder(spin_system,@sp_acquire,parameters,'labframe');
    
    % Apodisation
    fid=apodisation(spin_system,fid,{{'exp',10}});
    
    % Fourier transform and addition
    spectrum=spectrum+W(n)*fftshift(fft(fid,parameters.zerofill));
    
    % Plotting
    plot_1d(spin_system,real(spectrum),parameters); drawnow();
    
end

end

