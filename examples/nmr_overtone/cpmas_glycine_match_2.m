% Cross-polarization experiment between protons and 14N overtone
% transition in glycine under MAS. Glycine quadrupolar tensor da-
% ta comes from the paper by O'Dell and Ratcliffe:
%
%         http://dx.doi.org/10.1016/j.cplett.2011.08.030
%
% Hartmann-Hahn condition profile with a rough powder grid, as
% a function of spinning rate and 1H RF power.
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function cpmas_glycine_match_2()

% System specification
sys.magnet=14.10220742; sys.isotopes={'14N','1H'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);
inter.coupling.matrix{2,2}=[];
inter.zeeman.scalar={32.4  0.0};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.00]};
                  
% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=1000;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Spectrum setup
parameters.max_rank=5;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='rep_2ang_200pts_oct';
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
parameters.rf_dur=1e-4;

% Proton RF power array
rf_powers=linspace(10e3,200e3,50);

% Spinning rate array
spin_rates=linspace(20e3,90e3,50);

% Matrix preallocation
intensities=zeros(numel(rf_powers),numel(spin_rates));

% Loop over powers and spinning rates
parfor nk=1:numel(rf_powers)*numel(spin_rates)
    
    % Get the indices
    [n,k]=ind2sub([numel(rf_powers) numel(spin_rates)],nk);
    
    % Get a copy of parameters
    localpar=parameters;
    
    % Power setting
    localpar.rf_pwr=2*pi*[55e3 rf_powers(n)]/sin(theta);
    
    % Spinning rate
    localpar.rate=-spin_rates(k);
    
    % RF frequency on nitrogen OT
    localpar.rf_frq=8e3-2*localpar.rate;
    
    % Sweep width
    localpar.sweep=[localpar.rf_frq-4e3 ...
                    localpar.rf_frq+4e3];
    
    % Simulation
    spectrum=singlerot(spin_system,@overtone_cp,localpar,'qnmr');
    
    % Intensity
    intensities(nk)=sum(real(spectrum));
    
end

% Plot as image
figure();
imagesc([min(rf_powers)  max(rf_powers)]/1000,...
        [min(spin_rates) max(spin_rates)]/1000,intensities);
kxlabel('1H nutation frequency, kHz');
kylabel('Sample spinning rate, kHz');
set(gca,'YDir','normal'); colorbar;

end

