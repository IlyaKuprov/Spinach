% Multiple-quantum NMR experiment for a propanol 
% spin system. 
%
% Calculation time: seconds
% 
% m.g.concilio@soton.ac.uk
% jean-nicolas.dumez@univ-nantes.fr

function mqs_propanol()

% Magnet field
sys.magnet=14.1;

% Chemical shifts
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'};
inter.zeeman.scalar={3.438,3.438,1.429,1.429,0.775,0.775,0.775};

% 2J couplings
inter.coupling.scalar{1,3}=7.5;
inter.coupling.scalar{1,4}=7.5;
inter.coupling.scalar{3,5}=7.0;
inter.coupling.scalar{3,6}=7.0;
inter.coupling.scalar{3,7}=7.0;

% 3J couplings
inter.coupling.scalar{1,5}=0.5;
inter.coupling.scalar{1,6}=0.5;
inter.coupling.scalar{1,7}=0.5;
inter.coupling.scalar{7,7}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_group={'S2','S2','S3'};
bas.sym_spins={[1 2],[3 4],[5 6 7]};

% Algorithmic options
sys.enable={'greedy'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial and detection states
parameters.rho0=state(spin_system,'Lz','1H');
parameters.coil=state(spin_system,'L+','1H');

% Sequence parameters
parameters.angle=pi/2;
parameters.offset=[1200 1200];
parameters.sweep=[6000 2700];
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H','1H'};
parameters.axis_units='kHz';

% Tau delays
tau_max=[0.0333 0.0710 0.5000];

% Coherence to select
parameters.mqorder=+3;

% Get a figure going
figure(); scale_figure([2.5 1.5]);

% Loop over tau delays
for n=1:numel(tau_max)
    
    % tau delay    
    parameters.delay=tau_max(n);
             
    % Run simulation        
    fid=liquid(spin_system,@mqs,parameters,'nmr');
        
    % Apodization      
    fid=apodization(fid,'sinbell-2d');
        
    % Fourier transform      
    spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                               parameters.zerofill(1)));
                        
    % Get subplots for differnt tau values         
    subplot(1,3,n); 
    plot_2d(spin_system,abs(spectrum),parameters,20,...
            [0.02 0.2 0.02 0.2],2,256,6,'both'); 
    title(['\tau = ' num2str(parameters.delay) ' s']);          
    drawnow(); 
    
end

end

