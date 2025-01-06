% Multiple-quantum (MQ) NMR experiment for a coupled system 
% of six spins.
%
% Calculation time: minutes
% 
% mariagrazia.concilio@sjtu.edu.cn
% jean-nicolas.dumez@univ-nantes.fr

function mqs_six_spin()

% Magnetic field
sys.magnet=14.1;

% Chemical shifts
sys.isotopes={'1H','1H','1H','1H','1H','1H'};
inter.zeeman.scalar={-1.0 -0.5 0.0 +0.3 +0.7 +1.1};

% 3J couplings  
inter.coupling.scalar{1,2}=8.77;
inter.coupling.scalar{2,3}=7.40;
inter.coupling.scalar{3,4}=7.40;
inter.coupling.scalar{4,5}=8.77;
inter.coupling.scalar{5,6}=8.77;

% 4J couplings  
inter.coupling.scalar{1,3}=2.60;
inter.coupling.scalar{2,4}=2.60;
inter.coupling.scalar{3,5}=2.60;
inter.coupling.scalar{1,6}=2.60;
inter.coupling.scalar{5,6}=2.60;

% 5J couplings 
inter.coupling.scalar{1,4}=2.00;
inter.coupling.scalar{2,5}=2.00;
inter.coupling.scalar{2,6}=2.00;
inter.coupling.scalar{4,6}=2.00;

% 6J couplings
inter.coupling.scalar{1,5}=2.00;
inter.coupling.scalar{3,6}=2.00;
inter.coupling.scalar{6,6}=0;

% Coherence to select
parameters.mqorder=[+6 -1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

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
parameters.offset=[0 0];       % Hz
parameters.sweep=[2400 2400];  % Hz
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H','1H'};
parameters.axis_units='ppm';

% Tau delays
parameters.delay_1=0.1; % s
parameters.delay_2=0.1; % s
    
% Run simulation        
fid=liquid(spin_system,@mqs_refocus,parameters,'nmr');       
   
% Apodisation      
fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}});
   
% Fourier transform         
spectrum_conv=fftshift(fft2(fid,parameters.zerofill(2),...
                                parameters.zerofill(1)));
                        
% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum_conv),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'both');           
kxlabel('1Q / ppm'); kylabel('6Q / ppm');
    
end

