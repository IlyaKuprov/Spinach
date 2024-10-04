% Water protons in the unit cell of monohydrocalcite, inc-
% luding slow isotropic rotational diffusion and MAS. Fur-
% ther details in:
%
%       https://doi.org/10.1038/s41467-023-44381-x
%
% Calculation time: minutes, seconds with a GPU.
%
% guinevere.mathies@uni-konstanz.de

function sle_nmr_dd_csa_mhc()

% 400 MHz NMR
sys.magnet=9.4;

% Read CASTEP file
props=c2spinach('mhc.magres');

% Drop C, O, and Ca atoms 
drop_mask=ismember(props.symbols,{'C','O','Ca'});
props.symbols(drop_mask)=[];
props.std_geom(drop_mask,:)=[];
props.cst(drop_mask)=[];

% Keep two protons
sys.isotopes{1}='1H';
sys.isotopes{2}='1H';

% Convert shielding tensors into shift using the
% parametrisation of Huang et al. ACIE 2021
inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1};
inter.zeeman.matrix{2}=29.25*eye(3)-props.cst{4};

% Get coordinates
inter.coordinates{1}=props.std_geom(1,:);
inter.coordinates{2}=props.std_geom(4,:);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% This needs a GPU
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.rate=10000;
parameters.axis=[1 1 1];
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.spins={'1H'};
parameters.decouple={};
parameters.sweep=120000;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.offset=0;

% Isotropic rotational diffusion correlation times
tau_c=1e-6*[0.10 1.00 10.0 100.0 1000.0];

% Wigner function ranks 
max_rank=[2 3 5 7 13];

% Get figure going
figure(); hold on; 

% Loop over tau_c
for m=1:numel(tau_c)
    
    % Coorelation time and rank
    parameters.tau_c=tau_c(m);
    parameters.max_rank=max_rank(m);
    
    % Time-domain acquisition
    fid=gridfree(spin_system,@acquire,parameters,'nmr');

    % Apodisation
    fid=apodisation(spin_system,fid,{{'exp',6}});

    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));

    % Plotting
    plot_1d(spin_system,real(spectrum),parameters); 
    ylim([-1.0 10.0]); drawnow;
    
end

% Legend and labels
legend('\tau_c = 10^{-7} s',...
       '\tau_c = 10^{-6} s',...
       '\tau_c = 10^{-5} s',...
       '\tau_c = 10^{-4} s',...
       '\tau_c = 10^{-3} s','Location','NorthEast')
kylabel('amplitude, a.u.');

end

