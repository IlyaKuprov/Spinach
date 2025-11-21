% 13C MAS spectrum of glycine powder (assuming decoupling of 1H),
% computed using the Fokker-Planck MAS formalism and a spherical
% grid. The field dependence of the line shape of 13CA due to the
% presence of the quadrupolar 14N nucleus is shown. The calcula-
% tion is performed in the rotating frame with respect to 13C and
% the laboratory frame with respect to 14N.
%
% Calculation time: seconds.
%
% sanjay.vinod-kumar@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function mas_powder_gly_13c()

% Read CASTEP file
props=c2spinach('glycine.magres');

% Drop H and O atoms
drop_mask=ismember(props.symbols,{'H','O'});
props.symbols(drop_mask)=[];
props.std_geom(drop_mask,:)=[];
props.cst(drop_mask)=[];
props.efg(drop_mask)=[];

% Keep 13CA and 14N
sys.isotopes{1}='13C';
sys.isotopes{2}='14N';

% Convert shielding tensors into shift 
inter.zeeman.matrix{1}=-props.cst{2};
inter.zeeman.matrix{2}=-props.cst{5};

% Set isotropic chemical shifts to experimental values
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2],[43.6 110.0]);

% Cartesian coordinates
inter.coordinates={props.std_geom(2,:);
                   props.std_geom(5,:)};

% Quadrupolar interaction
inter.coupling.matrix=cell(2,2);
nqi=castep2nqi(props.efg{5},20.44e-3,1);
inter.coupling.matrix{2,2}=remtrace(nqi);

% Experiment setup
parameters.rate=10000;
parameters.axis=[1 1 1];
parameters.max_rank=5;
parameters.npoints=128;
parameters.zerofill=1024;
parameters.spins={'13C'};
parameters.grid='rep_2ang_200pts_sph';
parameters.verbose=0;

% Numerical rotating frame transforms
parameters.rframes={{'13C',1},{'14N',2}};

% Vary the field
fields=[4.7 9.4 14.1];

% Get a figure going
kfigure(); hold on;

for m=1:numel(fields)
					
    % Magnet field
    sys.magnet=fields(m);    
    
    % Basis set
    bas.formalism='sphten-liouv';
    bas.approximation='none';
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Experiment setup
    parameters.sweep=200*sys.magnet;      % unchanged in ppm
    parameters.offset=43.6*sys.magnet*spin('13C')/(2*pi*1e6);
    parameters.rho0=state(spin_system,'L+','13C');
    parameters.coil=state(spin_system,'L+','13C');
        
    % Lab frame Hamiltonian, then numerical rotating frames
    fid=singlerot(spin_system,@acquire,parameters,'labframe');
    
    % Apodisation
    fid=apodisation(spin_system,fid,{{'exp',6}});
    
    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));
    
    % Plotting
    plot_1d(spin_system,real(spectrum),parameters);

end

% Legend
klegend({'4.7 Tesla','9.4 Tesla','14.1 Tesla'},...
        'Location','NorthEast');

end

