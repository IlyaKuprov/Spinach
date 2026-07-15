% Methyl trosy in a rapidly rotating 13CH3 group
% of a slowly tumbling protein, simulated using 
% the Fokker-Planck formalism.
%
% ilya.kuprov@weizmann.ac.il

function trosy_methyl()

% Magnet field
sys.magnet=14.1;

% Cartesian coordinates
me_xyz=[-1.290651 1.279824 -0.248354;  % C
        -1.237804 1.209518 -1.347802;  % H
        -2.348228 1.329796  0.051440;  % H
        -0.799560 2.212894  0.057684]; % H

% Absolute shielding tensors (DFT)
shielding=cell(1,4);
shielding{1}=[166.0   -3.0    8.5;
               -3.0  178.6   -5.6;
                8.5   -5.6  165.0];
shielding{2}=[ 27.8   -0.9    1.3;
               -0.9   29.1   -1.9;
                1.3   -1.9   34.5];
shielding{3}=[ 34.6   -2.6   -1.1;
               -2.6   29.3   -0.2;
               -1.1   -0.2   26.7];
shielding{4}=[ 28.8    1.3    1.0;
                1.3   34.9    1.4;
                1.0    1.4   25.7];

% Convert shielding tensors into 
% chemical shift tensors and put
% them on resonance
shift_tensor=cell(1,4);
for n=1:4
    shift_tensor{n}=-remtrace(shielding{n});
end

% Methyl proton chemical shifts (guess)
shift_tensor{2}=shift_tensor{2}+0.8*eye(3);
shift_tensor{3}=shift_tensor{3}+1.0*eye(3);
shift_tensor{4}=shift_tensor{4}+1.2*eye(3);

% J-couplings
j_coupling=cell(4,4);
j_coupling{1,2}=125;
j_coupling{1,3}=125;
j_coupling{1,4}=125;
j_coupling{2,3}=-12;
j_coupling{2,4}=-12;
j_coupling{3,4}=-12;

% Spin system instances for the three methyl rotamers
sys.isotopes=repmat({'13C','1H','1H','1H'},1,3);
inter.coordinates=cell(1,12);
inter.zeeman.matrix=cell(1,12);
inter.coupling.scalar=cell(12,12);

% Assign spins to rotamers
inter.chem.parts={[1 2 3 4], [5 6 7 8], [9 10 11 12]};

% Equal rotamer populations
inter.chem.concs=[1 1 1];

% Interactions within methyl rotamer A
inter.coordinates{1}=me_xyz(1,:);
inter.coordinates{2}=me_xyz(2,:);
inter.coordinates{3}=me_xyz(3,:);
inter.coordinates{4}=me_xyz(4,:);
inter.zeeman.matrix{1}=shift_tensor{1};
inter.zeeman.matrix{2}=shift_tensor{2};
inter.zeeman.matrix{3}=shift_tensor{3};
inter.zeeman.matrix{4}=shift_tensor{4};
inter.coupling.scalar(1:4,1:4)=j_coupling;

% Interactions within methyl rotamer B
inter.coordinates{5}=me_xyz(1,:);
inter.coordinates{6}=me_xyz(4,:);
inter.coordinates{7}=me_xyz(2,:);
inter.coordinates{8}=me_xyz(3,:);
inter.zeeman.matrix{5}=shift_tensor{1};
inter.zeeman.matrix{6}=shift_tensor{4};
inter.zeeman.matrix{7}=shift_tensor{2};
inter.zeeman.matrix{8}=shift_tensor{3};
inter.coupling.scalar(5:8,5:8)=j_coupling;

% Interactions within methyl rotamer C
inter.coordinates{9}=me_xyz(1,:);
inter.coordinates{10}=me_xyz(3,:);
inter.coordinates{11}=me_xyz(4,:);
inter.coordinates{12}=me_xyz(2,:);
inter.zeeman.matrix{9}=shift_tensor{1};
inter.zeeman.matrix{10}=shift_tensor{3};
inter.zeeman.matrix{11}=shift_tensor{4};
inter.zeeman.matrix{12}=shift_tensor{2};
inter.coupling.scalar(9:12,9:12)=j_coupling;

% Formalism and basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Methyl turning generator
tau_m=1e-11; k_jump=1/(2*tau_m);
inter.chem.rates=k_jump*[-2  1  1;
                          1 -2  1;
                          1  1 -2];

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get a figure going
kfigure(); scale_figure([1.5 0.8]);

% Parameters: 13C NMR
parameters.tau_c=50e-9;
parameters.max_rank=3;
parameters.rho0=state(spin_system,'L+','13C');
parameters.coil=state(spin_system,'L+','13C');
parameters.decouple={};
parameters.spins={'13C'};
parameters.sweep=[-300 300];
parameters.npoints=1024;
parameters.axis_units='Hz';

% Call frequency domain detection
spectrum=gridfree(spin_system,@slowpass,parameters,'nmr');

% Do the plotting
subplot(1,2,1); plot_1d(spin_system,real(spectrum),parameters);

% Parameters: 1H NMR
parameters.tau_c=50e-9;
parameters.max_rank=3;
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.spins={'1H'};
parameters.sweep=[200 1000];
parameters.npoints=2048;
parameters.axis_units='Hz';

% Call frequency domain detection
spectrum=gridfree(spin_system,@slowpass,parameters,'nmr');

% Do the plotting
subplot(1,2,2); plot_1d(spin_system,real(spectrum),parameters);

end

