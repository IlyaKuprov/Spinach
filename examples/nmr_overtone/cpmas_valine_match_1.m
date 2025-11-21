% Cross-polarization experiment between protons and 14N overtone
% transition in N-acetylvaline under MAS, computed using Fokker-
% Planck formalism. Valine quadrupolar tensor data comes from:
%
%            http://dx.doi.org/10.1039/c4cp03994g
% 
% Hartmann-Hahn condition profile with a rough powder grid, as a
% function of 1H RF power.
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function cpmas_valine_match_1()

% System specification
sys.magnet=14.10220742; sys.isotopes={'14N','1H'};
inter.coupling.matrix{1,1}=eeqq2nqi(3.21e6,0.27,1,[0 0 0]);
inter.zeeman.eigs={[57.5 81.0 227.0],[0 0 0]};
inter.zeeman.euler={[-90 -90 -17],[0 0 0]};
inter.coupling.matrix{2,2}=[];
inter.coordinates={[0.00 0.00 0.00]
                   [1.00 0.00 0.00]};

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=2000;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Spectrum setup
parameters.max_rank=9;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=-19840;
parameters.grid='rep_2ang_800pts_sph';
parameters.sweep=[70.0e3 105.0e3];
parameters.npoints=256;
parameters.zerofill=256;
parameters.rho0=cos(theta)*state(spin_system,'Lz','1H')+...
                sin(theta)*state(spin_system,'Lx','1H');
parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+...
                sin(theta)*state(spin_system,'Lx','14N');
parameters.Nx=cos(theta)*operator(spin_system,'Lz','14N')+...
              sin(theta)*operator(spin_system,'Lx','14N');
parameters.Hx=cos(theta)*operator(spin_system,'Lz','1H')+...
              sin(theta)*operator(spin_system,'Lx','1H');
parameters.spins={'14N'};
parameters.method='average';
parameters.axis_units='kHz';
parameters.rf_frq=86.30e3;
parameters.rf_dur=7e-5;

% Proton RF power array
rf_powers=linspace(26e3,40e3,15);

% Start a new figure
kfigure(); scale_figure([2.5 1.0]);

% Proton RF power scan
for n=1:15
    
    % Subplot selection
    subplot(1,15,n);
    
    % Power setting
    parameters.rf_pwr=2*pi*[55e3 rf_powers(n)]/sin(theta);
    
    % Simulation
    spectrum=singlerot(spin_system,@overtone_cp,parameters,'qnmr');
       
    % Plotting
    plot_1d(spin_system,real(spectrum),parameters);
    axis([70 105 -2.5e-8 2.5e-8]); set(gca,'YTick',[]);
    ktitle([num2str(rf_powers(n)/1e3) ' kHz']); 
    kxlabel(''); set(gca,'XTick',[]); drawnow(); 
    
end

end

