% Fitting of 1H NMR spectrum of syn-2,4-difluoropentane with
% respect to J-couplings. See our paper for further details:
%
%        https://doi.org/doi/10.1021/acs.joc.4c00670
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% neil.wells@soton.ac.uk
% bruno.linclau@soton.ac.uk

function syn_difluoropentane()

% Load experimental data
load('syn_dfp_fluorine.mat','axis_ppm','spec'); axis_expt_f=axis_ppm;  spec_expt_f=spec;
load('syn_dfp_proton_a.mat','axis_ppm','spec'); axis_expt_ha=axis_ppm; spec_expt_ha=spec;
load('syn_dfp_proton_b.mat','axis_ppm','spec'); axis_expt_hb=axis_ppm; spec_expt_hb=spec;

% Normalise and shift the data
spec_expt_f=  7*spec_expt_f/max(spec_expt_f);
spec_expt_ha= 4*spec_expt_ha/max(spec_expt_ha)-0.1;
spec_expt_hb=10*spec_expt_hb/max(spec_expt_hb);

% Set the guess
guess=[ 6.2780   23.5067    5.1361    7.0537   24.9899   16.9264...
       48.0475    1.6017  -14.5413    0.7437    0.8771    0.6311];

% Set optimiser options
options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf);

% Get the figure going
figure(); scale_figure([2.00 1.75]);

% Run the optimisation
answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f,...
                             axis_expt_ha,spec_expt_ha,...
                             axis_expt_hb,spec_expt_hb,x),guess,options);

% Display the result
disp(answer);

end

% Least squares error function
function err=errfun(axis_expt_f, spec_expt_f,...
                    axis_expt_ha,spec_expt_ha,...
                    axis_expt_hb,spec_expt_hb,params)

% Silence Spinach
sys.output='hush';
sys.disable={'hygiene'};

% Run on GPU
sys.enable={'gpu'};

% Magnet induction 
sys.magnet=11.7464;

% Isotopes
sys.isotopes={'1H','1H','1H','1H','19F','1H','1H','1H','19F','1H','1H','1H'};  

% Chemical shifts
inter.zeeman.scalar{1}= 1.0189;
inter.zeeman.scalar{2}= 1.0189;
inter.zeeman.scalar{3}= 1.0189;
inter.zeeman.scalar{10}=1.0189;
inter.zeeman.scalar{11}=1.0189;
inter.zeeman.scalar{12}=1.0189;
inter.zeeman.scalar{4}=4.8515;
inter.zeeman.scalar{8}=4.8515;
inter.zeeman.scalar{5}=-173.4592;
inter.zeeman.scalar{9}=-173.4592;
inter.zeeman.scalar{6}=1.7910;
inter.zeeman.scalar{7}=2.1519;

% J-couplings
inter.coupling.scalar=cell(12,12);  
inter.coupling.scalar{1,4}=params(1);
inter.coupling.scalar{2,4}=params(1);
inter.coupling.scalar{3,4}=params(1);
inter.coupling.scalar{8,10}=params(1);
inter.coupling.scalar{8,11}=params(1);
inter.coupling.scalar{8,12}=params(1);
inter.coupling.scalar{1,5}=params(2);
inter.coupling.scalar{2,5}=params(2);
inter.coupling.scalar{3,5}=params(2);
inter.coupling.scalar{9,10}=params(2);
inter.coupling.scalar{9,11}=params(2);
inter.coupling.scalar{9,12}=params(2);
inter.coupling.scalar{4,6}=params(3);
inter.coupling.scalar{8,6}=params(3);
inter.coupling.scalar{4,7}=params(4);
inter.coupling.scalar{8,7}=params(4);
inter.coupling.scalar{5,6}=params(5);
inter.coupling.scalar{9,6}=params(5);
inter.coupling.scalar{5,7}=params(6);
inter.coupling.scalar{9,7}=params(6);
inter.coupling.scalar{4,5}=params(7);
inter.coupling.scalar{8,9}=params(7);
inter.coupling.scalar{5,9}=params(8);
inter.coupling.scalar{6,7}=params(9);

% Basis set and symmetry
bas.formalism='zeeman-hilb';
bas.approximation='none';
bas.sym_group={'S3','S3'};
bas.sym_spins={[1 2 3],[10 11 12]};

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas); 

% Sequence parameters (19F)
parameters_f.spins={'19F'};
parameters_f.rho0=state(spin_system,'L+','19F');
parameters_f.coil=state(spin_system,'L+','19F');
parameters_f.decouple={};
parameters_f.offset=-81655;
parameters_f.sweep=300;
parameters_f.npoints=512;
parameters_f.zerofill=2048;
parameters_f.axis_units='ppm';
parameters_f.invert_axis=1;

% Sequence parameters (1H A)
parameters_ha.spins={'1H'};
parameters_ha.rho0=state(spin_system,{'L+'},{4})+state(spin_system,{'L+'},{8});
parameters_ha.coil=state(spin_system,{'L+'},{4})+state(spin_system,{'L+'},{8});
parameters_ha.decouple={};
parameters_ha.offset=2426;
parameters_ha.sweep=128;
parameters_ha.npoints=256;
parameters_ha.zerofill=1024;
parameters_ha.axis_units='ppm';
parameters_ha.invert_axis=1;

% Sequence parameters (1H B)
parameters_hb.spins={'1H'};
parameters_hb.rho0=state(spin_system,{'L+'},{6})+state(spin_system,{'L+'},{7});
parameters_hb.coil=state(spin_system,{'L+'},{6})+state(spin_system,{'L+'},{7});
parameters_hb.decouple={};
parameters_hb.offset=986;
parameters_hb.sweep=350;
parameters_hb.npoints=512;
parameters_hb.zerofill=2048;
parameters_hb.axis_units='ppm';
parameters_hb.invert_axis=1;

% Simulation
fid_f=liquid(spin_system,@acquire,parameters_f,'nmr');
fid_ha=liquid(spin_system,@acquire,parameters_ha,'nmr');
fid_hb=liquid(spin_system,@acquire,parameters_hb,'nmr');

% Apodisation and scaling
fid_f= params(10)*apodization(fid_f, 'gaussian-1d',7.0)/4e3;
fid_ha=params(11)*apodization(fid_ha,'gaussian-1d',7.0)/4e3;
fid_hb=params(12)*apodization(fid_hb,'gaussian-1d',6.0)/4e3;

% Fourier transform
spec_theo_f= real(fftshift(fft(fid_f, parameters_f.zerofill)));  spec_theo_f= spec_theo_f(end:-1:1);
spec_theo_ha=real(fftshift(fft(fid_ha,parameters_ha.zerofill))); spec_theo_ha=spec_theo_ha(end:-1:1);
spec_theo_hb=real(fftshift(fft(fid_hb,parameters_hb.zerofill))); spec_theo_hb=spec_theo_hb(end:-1:1);

% Axis unification and interpolation
axis_theo_f= sweep2ticks(parameters_f.offset, parameters_f.sweep, parameters_f.zerofill);
axis_theo_ha=sweep2ticks(parameters_ha.offset,parameters_ha.sweep,parameters_ha.zerofill);
axis_theo_hb=sweep2ticks(parameters_hb.offset,parameters_hb.sweep,parameters_hb.zerofill);
axis_theo_f= -2*pi*1e6*axis_theo_f/spin_system.inter.basefrqs(5);
axis_theo_ha=-2*pi*1e6*axis_theo_ha/spin_system.inter.basefrqs(1);
axis_theo_hb=-2*pi*1e6*axis_theo_hb/spin_system.inter.basefrqs(1);
spec_theo_f= interp1(axis_theo_f,spec_theo_f,axis_expt_f,'pchip');    axis_theo_f=axis_expt_f;
spec_theo_ha=interp1(axis_theo_ha,spec_theo_ha,axis_expt_ha,'pchip'); axis_theo_ha=axis_expt_ha;
spec_theo_hb=interp1(axis_theo_hb,spec_theo_hb,axis_expt_hb,'pchip'); axis_theo_hb=axis_expt_hb;

% Plotting
subplot(3,1,1); plot(axis_expt_f,spec_expt_f,'ro','MarkerSize',1); 
hold on; plot(axis_theo_f,spec_theo_f,'b-'); hold off; axis tight; 
set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm'); kgrid;
subplot(3,1,2); plot(axis_expt_ha,spec_expt_ha,'ro','MarkerSize',1); 
hold on; plot(axis_theo_ha,spec_theo_ha,'b-'); hold off; axis tight; 
set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm'); kgrid;
subplot(3,1,3); plot(axis_expt_hb,spec_expt_hb,'ro','MarkerSize',1); 
hold on; plot(axis_theo_hb,spec_theo_hb,'b-'); hold off; axis tight; 
set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm'); kgrid;
klegend({'experiment','simulation'}); disp(params); drawnow;

% Error functional
err=norm(spec_expt_f-spec_theo_f)^2+...
    norm(spec_expt_ha-spec_theo_ha)^2+...
    norm(spec_expt_hb-spec_theo_hb)^2;

end

