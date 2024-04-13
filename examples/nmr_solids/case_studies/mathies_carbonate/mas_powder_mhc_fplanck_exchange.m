% Water protons in the unit cell of monohydrocalcite, 
% including position exchange and MAS.
%
% Calculation time: seconds.
%
% guinevere.mathies@uni-konstanz.de

function mas_powder_mhc_fplanck_exchange()

% 400 MHz NMR
sys.magnet=9.4;

% Read CASTEP file
props=c2spinach('mhc.magres');

% Drop C, O, and Ca atoms 
drop_mask=ismember(props.symbols,{'C','O','Ca'});
props.symbols(drop_mask)=[];
props.std_geom(drop_mask,:)=[];
props.cst(drop_mask)=[];

% Two reaction endpoints with two protons 
% each, swapped by the reaction
sys.isotopes{1}='1H';
sys.isotopes{2}='1H';
sys.isotopes{3}='1H';
sys.isotopes{4}='1H';

% Convert shielding tensors into shift using the
% parametrisation of Huang et al. ACIE 2021
inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1};
inter.zeeman.matrix{2}=29.25*eye(3)-props.cst{4};
inter.zeeman.matrix{3}=29.25*eye(3)-props.cst{4};
inter.zeeman.matrix{4}=29.25*eye(3)-props.cst{1};

% Get coordinates
inter.coordinates{1}=props.std_geom(1,:);
inter.coordinates{2}=props.std_geom(4,:);
inter.coordinates{3}=props.std_geom(4,:);
inter.coordinates{4}=props.std_geom(1,:);

% Chemical kinetics endpoints
inter.chem.parts={[1 2],[3 4]};

% Reaction rate matrix, Hz
inter.chem.rates=2e3*[-1  1; 
                       1 -1];

% Initial concentrations (arb. units)
inter.chem.concs=[1 1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Pulse-acquire setup
parameters.rate=10000;
parameters.axis=[1 1 1];
parameters.max_rank=13;
parameters.grid='rep_2ang_100pts_sph';
parameters.sweep=1/(5e-6);
parameters.npoints=512;
parameters.zerofill=1024;
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');

% Simulation
fid=singlerot(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

