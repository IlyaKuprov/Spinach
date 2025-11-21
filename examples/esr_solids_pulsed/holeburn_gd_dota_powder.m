% A hole burning simulation for a gadolinium ion. The soft pulse is 
% simulated using Fokker-Planck formalism. Zero-field splitting dis-
% ribution is sampled using the statistical parameters reported in
% Figure 5 of Raitsimring et al, App. Mag. Res. 28, 281-295 (2005).
% A numerical powder grid and numerical second-order rotating frame
% transformation are used.
%
% Note: non-central transition Gd(III) holes are very shallow.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function holeburn_gd_dota_powder()

% Initialize the spectra
spectrum_a=zeros(2048,1,'like',1i);
spectrum_b=zeros(2048,1,'like',1i);

% Get the sampling
[D,E,W]=zfs_sampling(30,5,1e-2); drawnow();

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
    parameters.rframes={{'E8',2}};
    parameters.verbose=0;
    
    % Soft pulse parameters
    parameters.pulse_rnk=2;
    parameters.pulse_phi=-pi/2;
    parameters.pulse_frq=-0.5e9;
    parameters.pulse_dur=50.0e-9;
    parameters.method='expm';
    
    % Simulation A
    parameters.pulse_pwr=0;
    fid_a=powder(spin_system,@holeburn,parameters,'labframe');
    
    % Simulation B
    parameters.pulse_pwr=2*pi*1e7;
    fid_b=powder(spin_system,@holeburn,parameters,'labframe');
    
    % Apodisation
    fid_a=apodisation(spin_system,fid_a,{{'exp',10}});
    fid_b=apodisation(spin_system,fid_b,{{'exp',10}});
    
    % Fourier transform and addition
    spectrum_a=spectrum_a+W(n)*fftshift(fft(fid_a,parameters.zerofill));
    spectrum_b=spectrum_b+W(n)*fftshift(fft(fid_b,parameters.zerofill));
    
    % Plotting
    plot_1d(spin_system,real(spectrum_a),parameters,'r-'); hold on;
    plot_1d(spin_system,real(spectrum_b),parameters,'b-'); hold off; drawnow;
    
end

end

