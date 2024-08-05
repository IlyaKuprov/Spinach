% Fitting of a 79Br MAS NMR spectrum of potassium bromide 
% with respect to the quadrupole coupling constant.
%
% The spectrum cannot be fitted with a single quadrupolar
% tensor; at least 3 are necessary, likely due to a dist-
% ribution of electrostatic environments in the powder.
%
% Calculation time: hours.
%
% sanjay.vinod-kumar@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% i.kuprov@soton.ac.uk

function kbr_mas_fitting()

% Load and normalise the data
kbr_data=readmatrix('KBr_400MHz_2kHz.txt');
spec_expt=kbr_data(:,2)/100;

% Set instrumental variables
sys.magnet=9.3659;          % magnet field, Tesla
parameters.rate=2000;       % MAS rate, Hz
parameters.offset=6034.96;  % receiver offset, Hz
parameters.sweep=1e5;       % sweep width, Hz
parameters.npoints=4096;    % points to acquire
parameters.zerofill=32768;  % zerofill to

guess=[60.0933        ...  % chemical shift, ppm
       13.7569 1.6424 ...  % three sets of XX, YY NQI eigenvalues, kHz
       4.0779  4.5179 ...
       1.5885  0.9449 ...
       263.9835       ...  % relaxation rate (Hz)
       40 32 28 ];         % component weights

% Set optimizer options
options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf);

% Get a figure going
figure(); scale_figure([1.75 1.50]);

% Run the optimisation
fminsearch(@errfun,guess,options);

    % Least squares error function
    function err=errfun(params)

        % Silence Spinach
        sys.output='hush';
        sys.disable={'hygiene','trajlevel'};
        
        % Spin system
        sys.isotopes={'79Br','79Br','79Br'};

        % Chemical shift, ppm
        inter.zeeman.scalar={params(1) params(1) params(1)};

        % Quadrupolar coupling
        inter.coupling.matrix{1,1}=1e3*diag([params(2) params(3) -(params(2)+params(3))]);
        inter.coupling.matrix{2,2}=1e3*diag([params(4) params(5) -(params(4)+params(5))]);
        inter.coupling.matrix{3,3}=1e3*diag([params(6) params(7) -(params(6)+params(7))]);

        % Basis set
        bas.formalism='sphten-liouv';
        bas.approximation='IK-0';
        bas.level=1; bas.projections=+1;

        % Relaxation theory
        inter.relaxation={'t1_t2'};
        inter.r1_rates={params(8) params(8) params(8)};
        inter.r2_rates={params(8) params(8) params(8)};
        inter.equilibrium='zero';
        inter.rlx_keep='diagonal';

        % Spinach housekeeping
        spin_system=create(sys,inter);
        spin_system=basis(spin_system,bas);

        % Experiment parameters
        parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
        parameters.max_rank=50;
        parameters.grid='rep_2ang_200pts_oct';
        parameters.spins={'79Br'};
        parameters.decouple={};
        parameters.axis_units='Hz';
        parameters.invert_axis=1;
        parameters.rho0=params(9) *state(spin_system,{'L+'},{1})+...
                        params(10)*state(spin_system,{'L+'},{2})+...
                        params(11)*state(spin_system,{'L+'},{3});
        parameters.coil=state(spin_system,'L+','79Br');

        % Simulation
        fid=singlerot(spin_system,@acquire,parameters,'nmr');

        % Apodisation
        fid=apodization(fid,'none-1d');

        % Fourier transform and scaling
        spec_theo=real(fftshift(fft(fid,parameters.zerofill)));

        % Plotting
        plot_1d(spin_system,spec_expt,parameters,'r.'); hold on;
        plot_1d(spin_system,spec_theo,parameters,'b-'); hold off; 
        
        % Display the parameters
        disp(params); drawnow;

        % Error functional
        err=norm(spec_expt-spec_theo,2)^2;

    end

end

