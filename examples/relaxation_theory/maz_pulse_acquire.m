% Methylaziridine pulse-acquire, showing the effect of the scalar
% relaxation of the second kind due to the fast quadrupolar rela-
% xation of the 14N nuclei.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% tim.claridge@chem.ox.ac.uk
% barbara.odell@chem.ox.ac.uk

function maz_pulse_acquire()

% Magnet induction
sys.magnet=11.75;

% Isotopes
sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H'};

% Absolute shielding (vacuum DFT)
inter.zeeman.matrix{1}=[32.7075    6.1977   -3.6061
                         3.3757   29.2462   -0.8672
                        -6.8603   -3.1603   29.2097];    % H3
inter.zeeman.matrix{2}=[27.5505   -0.4604   -3.9670
                         2.4491   30.8072    3.5941
                         0.0951    5.9399   33.5623];    % H2
inter.zeeman.matrix{3}=[27.8325   -1.3662   -3.4664
                        -0.3195   23.8084   -0.9356
                         1.9567   -0.5013   39.0966];    % H4
inter.zeeman.matrix{4}=[257.9722    16.9286   -28.3305
                         22.7951   157.2701    18.5357
                        -13.8903   -20.6605   275.1264]; % N
inter.zeeman.matrix{5}=[28.3494    3.7981   -3.9745
                        -1.0065   31.8697   -4.2119
                        -1.1997  -10.1395   38.9642];    % HN
inter.zeeman.matrix{6}=[33.2953    4.3402    0.8099
                         3.3922   29.9675    1.2225
                        -0.7891    0.1926   26.6461];    % Me
inter.zeeman.matrix{7}=[33.2302   -4.7306    1.2366
                        -3.6419   30.0443   -1.1180
                         -0.3469   -0.7420   26.5551];   % Me
inter.zeeman.matrix{8}=[31.4875    0.1078   -2.9339
                         0.1534   27.9566    0.6704
                        -3.2188    0.7748   35.7251];    % Me

% Assign isotropic components from the experiment
expt_shift=[1.3 1.7 1.9 0.0 0.1 1.2 1.2 1.2];
for n=1:numel(inter.zeeman.matrix)
    inter.zeeman.matrix{n}=inter.zeeman.matrix{n}-eye(3)*trace(inter.zeeman.matrix{n})/3;
    inter.zeeman.matrix{n}=eye(3)*expt_shift(n)-inter.zeeman.matrix{n};
end
                        
% Quadrupole couplings (vacuum DFT)
inter.coupling.matrix=cell(8);
inter.coupling.matrix{4,4}=[-1.2932   0.6251   1.8700
                             0.6251   1.7170  -2.3127
                             1.8700  -2.3127  -0.4238]*1e6;
                         
% Scalar couplings (vacuum DFT)
inter.coupling.scalar=cell(8);
inter.coupling.scalar{1,2}=-11.6;
inter.coupling.scalar{1,3}=5.5;
inter.coupling.scalar{1,5}=7.5;
inter.coupling.scalar{2,3}=1.5;
inter.coupling.scalar{2,5}=12.7;
inter.coupling.scalar{3,6}=4.5;
inter.coupling.scalar{3,7}=4.5;
inter.coupling.scalar{3,8}=4.5;
inter.coupling.scalar{4,1}=4.4;
inter.coupling.scalar{4,3}=5.2;
inter.coupling.scalar{4,5}=44.8;
inter.coupling.scalar{5,6}=1.0;
inter.coupling.scalar{5,7}=1.0;
inter.coupling.scalar{5,8}=1.0;

% Coordinates (Angstrom, vacuum DFT)
inter.coordinates={[-1.812977   -1.098554    0.444452]
                   [-0.842255   -1.116468   -1.122812]
                   [ 0.176983   -0.023347    1.610045]
                   [-0.907675    0.811819   -0.038237]
                   [-0.590140    1.188036   -0.936471]
                   [ 2.093648    0.832363    0.088912]
                   [ 2.084298   -0.944150    0.159225]
                   [ 1.387021   -0.108214   -1.250363]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=3;

% Relaxation superoperator
inter.relaxation={'redfield','SRSK'};
inter.equilibrium='zero';
inter.srsk_sources=4;
inter.rlx_keep='secular';
inter.tau_c={200e-12};

% Algorithmic options
sys.tols.inter_cutoff=2.0;
sys.tols.prox_cutoff=4.0;
sys.disable={'krylov'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=500;
parameters.sweep=1400;
parameters.npoints=4096;
parameters.zerofill=16536;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

