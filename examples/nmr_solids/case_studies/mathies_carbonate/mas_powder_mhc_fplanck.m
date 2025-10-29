% All protons in the unit cell of monohydrocalcite, magic
% angle spinning NMR simulation. Further details in:
%
%      https://doi.org/10.1038/s41467-023-44381-x
%
% Calculation time: hours, minutes on a GPU.
%
% guinevere.mathies@uni-konstanz.de

function mas_powder_mhc_fplanck()

% 400 MHz NMR
sys.magnet=9.4;

% Read CASTEP file
props=c2spinach('mhc.magres');

% Drop C, O, and Ca atoms 
drop_mask=ismember(props.symbols,{'C','O','Ca'});
props.symbols(drop_mask)=[];
props.std_geom(drop_mask,:)=[];
props.cst(drop_mask)=[];

% Get isotopes
sys.isotopes={};
for n=1:numel(props.symbols)
    if strcmp(props.symbols{n},'H')
        sys.isotopes{n}='1H';
    else
        error('unexpected atom.');
    end
end

% Convert shielding tensors into shift using the
% parametrisation of Huang et al. ACIE 2021
inter.zeeman.matrix=cell(1,numel(props.cst));
for n=1:numel(props.cst)
    if strcmp(sys.isotopes{n},'1H')
        inter.zeeman.matrix{n}=29.25*eye(3)-props.cst{n};
    end
end

% Get coordinates
inter.coordinates=mat2cell(props.std_geom,ones(18,1));

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=3; 
bas.projections=+1;

% Interaction cut-off, Hz
sys.tols.inter_cutoff=500;

% This needs a GPU
% sys.enable={'gpu'};

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
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.verbose=1;

% Simulation
fid=singlerot(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

