% Methyl trosy in a rapidly rotating 13CH3 group
% of a slowly tumbling protein, simulated using 
% the Fokker-Planck formalism. Alanine is used
% as a geometric model.
%
% ilya.kuprov@weizmann.ac.il

function trosy_methyl()

% Magnet field
sys_big.magnet=14.1;

% Cartesian coordinates: CA carbon
ca_xyz=[-0.584588 0.114350 0.394793];

% Cartesian coordinates: Me group
me_xyz=[-1.290651 1.279824 -0.248354;
        -1.237804 1.209518 -1.347802;
        -2.348228 1.329796 0.051440;
        -0.799560 2.212894 0.057684];

% Cartesian coordinates: Me carbon
cb_xyz=me_xyz(1,:);

% Methyl group turning axis
turning_axis=cb_xyz-ca_xyz;

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
% chemical shift tensors
shift_ref=[188.0 31.8 31.8 31.8];
inter.coordinates=num2cell(me_xyz,2);
inter.zeeman.matrix=cell(1,5);
for n=1:4

    % Convert chemical shielding to chemcial shift
    inter.zeeman.matrix{n}=shift_ref(n)*eye(3)-shielding{n};

    % Put chemical shift on resonance with transmitter
    inter.zeeman.matrix{n}=remtrace(inter.zeeman.matrix{n});

end

% J-couplings for the methyl group
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2}=125;
inter.coupling.scalar{1,3}=125;
inter.coupling.scalar{1,4}=125;
inter.coupling.scalar{2,3}=-12;
inter.coupling.scalar{2,4}=-12;
inter.coupling.scalar{3,4}=-12;

% Spin system instances for the three methyl rotamers
sys_big.isotopes=repmat({'13C','1H','1H','1H'},1,3);
inter_big.coordinates=cell(1,12);
inter_big.zeeman.matrix=cell(1,12);
inter_big.coupling.scalar=cell(12,12);

% Rotamers
for r=1:3

    % Get the methyl group turning matrix
    R_me=anax2dcm(turning_axis,2*pi*(r-1)/3);

    % Turn the methyl group
    current_xyz=(me_xyz-cb_xyz)*R_me'+cb_xyz;

    % Assign atoms
    for n=1:4
        inter_big.coordinates{(r-1)*4+n}=current_xyz(n,:);
    end

    % Assign CSAs
    for n=1:4
        inter_big.zeeman.matrix{(r-1)*4+n}=R_me*inter.zeeman.matrix{n}*R_me';
    end

    % Assign J-couplings
    for n=1:4
        for m=1:4
            inter_big.coupling.scalar{(r-1)*4+n,(r-1)*4+m}=inter.coupling.scalar{n,m};
        end
    end

end

% Formalism and basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
bas.zero_quantum={'1H'};

% Assign spins to rotamers
inter_big.chem.parts=cell(1,3);
for n=1:3
    inter_big.chem.parts{n}=4*(n-1)+(1:4);
end

% Equal rotamer populations
inter_big.chem.concs=[1 1 1];

% Methyl jump generator
tau_m=1e-11; k_jump=1/(2*tau_m);
inter_big.chem.rates=k_jump*[-2  1  1;
                              1 -2  1;
                              1  1 -2];

% Spinach housekeeping
spin_system=create(sys_big,inter_big);
spin_system=basis(spin_system,bas);

% Spectrum parameters
parameters.tau_c=50e-9;
parameters.max_rank=3;
parameters.rho0=state(spin_system,'L+','13C');
parameters.coil=state(spin_system,'L+','13C');
parameters.decouple={};
parameters.spins={'13C'};
parameters.sweep=[-300 300];
parameters.npoints=2048;
parameters.axis_units='Hz';

% Call frequency domain detection
spectrum=gridfree(spin_system,@slowpass,parameters,'nmr');

% Do the plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

