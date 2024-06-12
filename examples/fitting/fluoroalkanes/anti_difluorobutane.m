% Fitting of 1H NMR spectrum of anti-2,3-difluorobutane with
% respect to J-couplings. See our paper for further details:
%
%        https://doi.org/doi/10.1021/acs.joc.4c00670
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% neil.wells@soton.ac.uk
% bruno.linclau@soton.ac.uk

function anti_difluorobutane()

% Load experimental data
load('anti_dfb_proton.mat','ch_axis_hz','ch_expt_data',...
                           'me_axis_hz','me_expt_data');  

% Normalise the data
ch_expt_data=-2*ch_expt_data/trapz(ch_axis_hz,ch_expt_data); 
me_expt_data=-6*me_expt_data/trapz(me_axis_hz,me_expt_data);

% Concatenate spectral intervals
expt_data=[ch_expt_data; me_expt_data];
axis_hz=[ch_axis_hz; me_axis_hz];

% Broaden out Z1 shim problem
filter=exp(-10*linspace(-1,1,100).^2);
expt_data=conv(expt_data,filter,'same')/sum(filter);

% Set the guess
guess=[24.08  6.49  1.44  3.59  15.73  47.76  -13.58  26.56  1.7];

% Get a figure going
figure(); scale_figure([1.0 1.5]);

% Set optimiser options
options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf);

% Run the optimisation
answer=fminsearch(@(x)errfun(axis_hz,expt_data,x),guess,options);

% Display the result
disp(answer);

end

% Least squares error function
function err=errfun(axis_hz,expt_data,params)

% Silence Spinach
sys.output='hush';
sys.disable={'hygiene'};

% Absorb parameters
j_f_ch3_near=params(1);
j_h_ch3=params(2);
j_f_ch3_far=params(3);
j3_h_h=params(4);
j3_f_h=params(5);
j2_f_h=params(6);
j3_f_f=params(7);
lw=params(8); 
A=params(9);

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','19F','19F'};

% Magnet field
sys.magnet=11.7464;

% Chemical shifts
inter.zeeman.scalar{1}=1.32375;
inter.zeeman.scalar{2}=1.32375;
inter.zeeman.scalar{3}=1.32375;
inter.zeeman.scalar{4}=1.32375;
inter.zeeman.scalar{5}=1.32375;
inter.zeeman.scalar{6}=1.32375;
inter.zeeman.scalar{7}=4.61040;
inter.zeeman.scalar{8}=4.61040;
inter.zeeman.scalar{9}= 0.00;
inter.zeeman.scalar{10}=0.00;

% Scalar couplings
inter.coupling.scalar=cell(10,10);
inter.coupling.scalar{1,9}=j_f_ch3_near;
inter.coupling.scalar{2,9}=j_f_ch3_near;
inter.coupling.scalar{3,9}=j_f_ch3_near;
inter.coupling.scalar{4,10}=j_f_ch3_near;
inter.coupling.scalar{5,10}=j_f_ch3_near;
inter.coupling.scalar{6,10}=j_f_ch3_near;
inter.coupling.scalar{1,10}=j_f_ch3_far;
inter.coupling.scalar{2,10}=j_f_ch3_far;
inter.coupling.scalar{3,10}=j_f_ch3_far;
inter.coupling.scalar{4,9}=j_f_ch3_far;
inter.coupling.scalar{5,9}=j_f_ch3_far;
inter.coupling.scalar{6,9}=j_f_ch3_far;
inter.coupling.scalar{1,7}=j_h_ch3;
inter.coupling.scalar{2,7}=j_h_ch3;
inter.coupling.scalar{3,7}=j_h_ch3;
inter.coupling.scalar{4,8}=j_h_ch3;
inter.coupling.scalar{5,8}=j_h_ch3;
inter.coupling.scalar{6,8}=j_h_ch3;
inter.coupling.scalar{7,8}=j3_h_h;
inter.coupling.scalar{8,9}=j3_f_h;
inter.coupling.scalar{7,10}=j3_f_h;
inter.coupling.scalar{7,9}=j2_f_h;
inter.coupling.scalar{8,10}=j2_f_h;
inter.coupling.scalar{9,10}=j3_f_f;

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
parameters.offset=1500;
parameters.sweep=1800;
parameters.npoints=4096;
parameters.zerofill=32768;
parameters.axis_units='Hz';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'gaussian-1d',lw);

% Fourier transform
sim_spec=A*fftshift(fft(fid,parameters.zerofill))/1e6;
sim_spec=sim_spec(end:-1:1);

% Axis unification
sim_axis=sweep2ticks(parameters.offset,parameters.sweep,parameters.zerofill);
sim_spec=interp1(sim_axis,sim_spec,axis_hz,'pchip');

% Plotting
subplot(2,1,1); plot(axis_hz,real(expt_data),'ro','MarkerSize',1); hold on;
plot(axis_hz,real(sim_spec),'b-'); hold off; xlim([2255 2355]); kgrid;
kxlabel('Frequency, Hz'); klegend({'experiment','simulation'});
subplot(2,1,2); plot(axis_hz,real(expt_data),'ro','MarkerSize',1); hold on;
plot(axis_hz,real(sim_spec),'b-'); hold off; xlim([635 690]); kgrid;
kxlabel('Frequency, Hz'); klegend({'experiment','simulation'}); drawnow;

% Error functional
err=norm(real(expt_data)-real(sim_spec))^2;

end

