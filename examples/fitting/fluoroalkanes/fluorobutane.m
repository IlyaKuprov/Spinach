% Fitting of 1H NMR spectrum of 2-fluoropentane with respect
% to J-couplings. See our paper for further details:
%
%        https://doi.org/doi/10.1021/acs.joc.4c00670
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% neil.wells@soton.ac.uk
% bruno.linclau@soton.ac.uk

function fluorobutane()

% Load experimental data
load('fluorobutane_fluorine.mat','spec_f','axis_f'); 
load('fluorobutane_proton.mat','spec_ch', 'axis_ch',...
                               'spec_ch2','axis_ch2');

% Normalise the data
spec_ch= -1*spec_ch/trapz(axis_ch,spec_ch); 
spec_ch2=-2*spec_ch2/trapz(axis_ch2,spec_ch2);
spec_f=  -1*spec_f/trapz(axis_f,spec_f);

% Concatenate spectral intervals
spec_h=[spec_ch; spec_ch2];
axis_h=[axis_ch; axis_ch2];

% Set the guess
guess=[23.9529    6.2219    7.4903    7.1692    4.9608   17.4939...
       26.2802   48.6869  -14.0924    1.8099    4.3571    ];

% Set optimiser options
options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf);

% Get the figure going
figure(); scale_figure([2.00 1.75]);

% Run the optimisation
answer=fminsearch(@(x)errfun(axis_h,spec_h,axis_f,spec_f,x),guess,options);

% Display the result
disp(answer);

end

% Least squares error function
function err=errfun(axis_h,expt_h,axis_f,expt_f,params)

% Silence Spinach
sys.output='hush';
sys.disable={'hygiene'};

% Absorb parameters
j_f_ch3=params(1);
j_h_ch3_near=params(2);
j_h_ch3_far=params(3);
j3_h_h_one=params(4);
j3_h_h_two=params(5);
j3_f_h_one=params(6);
j3_f_h_two=params(7);
j2_f_h=params(8);
j2_h_h=params(9);
A=params(10);
B=params(11);

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','19F'};

% Magnet field
sys.magnet=11.7464;

% Chemical shifts
inter.zeeman.scalar{1}=0.982;
inter.zeeman.scalar{2}=0.982;     % CH3 far
inter.zeeman.scalar{3}=0.982;
inter.zeeman.scalar{4}=1.333;
inter.zeeman.scalar{5}=1.333;     % CH3 near
inter.zeeman.scalar{6}=1.333;
inter.zeeman.scalar{7}=4.6075;    % CHF
inter.zeeman.scalar{8}=1.5949;    % CH2 
inter.zeeman.scalar{9}=1.6890;    % CH2
inter.zeeman.scalar{10}=-173.184; % F

% Scalar couplings
inter.coupling.scalar=cell(10,10);
inter.coupling.scalar{1,9}=j_h_ch3_far;
inter.coupling.scalar{2,9}=j_h_ch3_far;
inter.coupling.scalar{3,9}=j_h_ch3_far;
inter.coupling.scalar{1,8}=j_h_ch3_far;
inter.coupling.scalar{2,8}=j_h_ch3_far;
inter.coupling.scalar{3,8}=j_h_ch3_far;
inter.coupling.scalar{4,7}=j_h_ch3_near;
inter.coupling.scalar{5,7}=j_h_ch3_near;
inter.coupling.scalar{6,7}=j_h_ch3_near;
inter.coupling.scalar{4,10}=j_f_ch3;
inter.coupling.scalar{5,10}=j_f_ch3;
inter.coupling.scalar{6,10}=j_f_ch3;
inter.coupling.scalar{7,8}=j3_h_h_two;
inter.coupling.scalar{7,9}=j3_h_h_one;
inter.coupling.scalar{7,10}=j2_f_h;
inter.coupling.scalar{8,9}=j2_h_h;
inter.coupling.scalar{8,10}=j3_f_h_two;
inter.coupling.scalar{9,10}=j3_f_h_one;

% Basis set and symmetry
bas.formalism='zeeman-hilb';
bas.approximation='none';
bas.sym_group={'S3','S3'};
bas.sym_spins={[1 2 3],[4 5 6]};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=1400;
parameters.sweep=2000;
parameters.npoints=4096;
parameters.zerofill=32768;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation and scaling
fid=apodization(fid,'gaussian-1d',6.0)/2e3;

% Fourier transform
sim_h=A*fftshift(fft(fid,parameters.zerofill));

% Axis generation
ax=linspace(-parameters.sweep/2,parameters.sweep/2,numel(sim_h))+parameters.offset;
ax=1000000*(2*pi)*ax/(spin(parameters.spins{1})*spin_system.inter.magnet);

% Axis unification
sim_h=interp1(ax,sim_h,axis_h,'pchip');

% Cut-out point
sim_h(3770:3820)=0;
expt_h(3770:3820)=0;

% Plotting
subplot(3,1,1); plot(axis_h,real(expt_h),'ro','MarkerSize',1); 
hold on; plot(axis_h,real(sim_h),'b-'); hold off; xlim([1.5 1.8]); 
kgrid; set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm');
subplot(3,1,2); plot(axis_h,real(expt_h),'ro','MarkerSize',1); 
hold on; plot(axis_h,real(sim_h),'b-'); hold off; xlim([4.51 4.70]);
kgrid; set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm');

% Sequence parameters
parameters.spins={'19F'};
parameters.rho0=state(spin_system,'L+','19F');
parameters.coil=state(spin_system,'L+','19F');
parameters.decouple={};
parameters.offset=-81529;
parameters.sweep=250;
parameters.npoints=512;
parameters.zerofill=2048;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation and scaling
fid=apodization(fid,'gaussian-1d',6.0)/1e3;

% Fourier transform
sim_f=B*fftshift(fft(fid,parameters.zerofill));

% Axis generation
ax=linspace(-parameters.sweep/2,parameters.sweep/2,numel(sim_f))+parameters.offset;
ax=1000000*(2*pi)*ax/(spin(parameters.spins{1})*spin_system.inter.magnet);

% Axis unification
sim_f=interp1(ax,sim_f,axis_f,'pchip');

% Plotting
subplot(3,1,3); plot(axis_f,real(expt_f),'r.','MarkerSize',1); 
hold on; plot(axis_f,real(sim_f),'b-'); hold off; xlim([-173.4 -172.95]); 
kgrid; set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm');

% Print the parameters
disp(params); drawnow;

% Error functional
err=norm(real(expt_h)-real(sim_h))^2+...
    norm(real(expt_f)-real(sim_f))^2;

end

