% Earth's field NMR Simulation for 2,6-difluoropyridine; replicates 
% simulated spectra in Figure 7 of 
%
%              https://doi.org/10.1016/j.jmr.2023.107540
%
% without the weighted addition of the uncoupled 1H signal.
%
% Calculation time: seconds.
%
% Adam Altenhof, adamaltenhof@gmail.com
% Derrick Kaseman, kaseman1@llnl.gov

function earth_field_dfp()

% Earth's field
sys.magnet=2*pi*2312.45/spin('1H');

% Isotopes and labels
sys.isotopes={'1H', '1H', '1H', '19F', '19F', '14N'};
sys.labels=  {'H3', 'H4', 'H5', 'F2',  'F6',  'N1' };

% Chemical shifts
inter.zeeman.scalar={6.98 8.06 6.98 -70.69 -70.69 0};

% J-couplings (experimental)
inter.coupling.scalar=cell(6,6);
inter.coupling.scalar{idxof(sys,'F2'),idxof(sys,'N1')} =  37.35;
inter.coupling.scalar{idxof(sys,'F6'),idxof(sys,'N1')} =  37.35;
inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'N1')} =  (1.03+0.71)/2; % Table 2 makes no sense
inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'N1')} =  (1.03+0.71)/2; % Table 2 makes no sense
inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'F2')} =  -2.47; 
inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'F6')} =  -2.47;
inter.coupling.scalar{idxof(sys,'H4'),idxof(sys,'F2')} =   8.08;
inter.coupling.scalar{idxof(sys,'H4'),idxof(sys,'F6')} =   8.08;
inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'F2')} =   1.29;
inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'F6')} =   1.29;
inter.coupling.scalar{idxof(sys,'F2'),idxof(sys,'F6')} = -12.23;
inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'H4')} =   7.92;
inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'H4')} =   7.92;
inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'H5')} =   0.55;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'t1_t2'};

inter.rlx_keep='diagonal';
inter.equilibrium='zero';

% Sequence parameters
parameters.sweep=1e4;
parameters.npoints=60000;
parameters.zerofill=2^19;
parameters.offset=0;
parameters.spins={'1H'};
parameters.axis_units='Hz';
parameters.invert_axis=0;
parameters.detection='uniaxial';
parameters.flip_angle=pi/2;

% This needs a GPU
sys.enable={'gpu'}; 

% 14N relaxation rates, Hz
R_14N=[0 10 100 500 1e3 1e4 1e5 1e6];

% Get the figure going
figure(); scale_figure([1.0 2.0]);

% Loop over 14N relaxation rates
for n=1:numel(R_14N)

    % Set relaxation rates
    inter.r1_rates={0.22 0.22 0.22 0.22 0.22 R_14N(n)};
    inter.r2_rates={0.22 0.22 0.22 0.22 0.22 R_14N(n)};

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
     
    % Simulation
    fid=liquid(spin_system,@zerofield,parameters,'labframe');
    
    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));
    spectrum=real(spectrum)/max(real(spectrum));

    % Frequency axis
    axis_hz=sweep2ticks(parameters.offset,...
                        parameters.sweep,...
                        parameters.zerofill);
    
    plot(axis_hz,spectrum-n/3); xlim([2303 2322]);
    kxlabel('Frequency, Hz'); set(gca,'YTick',{}); 
    kgrid; hold on; drawnow;

end

end

