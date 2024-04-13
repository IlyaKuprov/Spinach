% Cross-polarization experiment between protons and 14N overtone
% transition in N-acetylvaline under MAS, computed using Fokker-
% Planck formalism. Valine quadrupolar tensor data comes from:
%
%            http://dx.doi.org/10.1039/c4cp03994g
% 
% Hartmann-Hahn condition profile with a rough powder grid, as a
% function of 1H RF power and spinning rate.
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function cpmas_valine_match_2()

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
inter.damp_rate=10000;

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
                sin(theta)*(state(spin_system,'L+','1H')+...
                            state(spin_system,'L-','1H'))/2;
parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+...
                sin(theta)*(state(spin_system,'L+','14N')+...
                            state(spin_system,'L-','14N'))/2;
parameters.Nx=cos(theta)*operator(spin_system,'Lz','14N')+...
              sin(theta)*(operator(spin_system,'L+','14N')+...
                          operator(spin_system,'L-','14N'))/2;
parameters.Hx=cos(theta)*operator(spin_system,'Lz','1H')+...
              sin(theta)*(operator(spin_system,'L+','1H')+...
                          operator(spin_system,'L-','1H'))/2;
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
    localpar.rf_frq=46.30e3-2*localpar.rate;
    
    % Sweep width
    localpar.sweep=[localpar.rf_frq-15e3 ...
        localpar.rf_frq+15e3];
    
    % Simulation
    spectrum=singlerot(spin_system,@overtone_cp,localpar,'qnmr');
    
    % Intensity
    intensities(nk)=sum(real(spectrum));
    
end

% Plot as image
figure();
imagesc([min(rf_powers) max(rf_powers)],...
        [min(spin_rates) max(spin_rates)],intensities);
xlabel('1H RF power, Hz'); set(gca,'YDir','normal');
ylabel('Sample spinning rate, Hz');

end

