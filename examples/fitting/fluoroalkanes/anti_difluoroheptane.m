% Fitting of 1H NMR spectrum of anti-3,5-difluoroheptane with
% respect to J-couplings. See our paper for further details:
%
%        https://doi.org/doi/10.1021/acs.joc.4c00670
%
% Methyl groups are ghosted out because they do not influence
% the signals in question.
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% neil.wells@soton.ac.uk
% bruno.linclau@soton.ac.uk

function anti_difluoroheptane()

% Load experimental data
load('anti_dfh_fluorine.mat','axis_ppm','spec'); axis_expt_f =axis_ppm; spec_expt_f =spec;
load('anti_dfh_proton_a.mat','axis_ppm','spec'); axis_expt_ha=axis_ppm; spec_expt_ha=spec;
load('anti_dfh_proton_b.mat','axis_ppm','spec'); axis_expt_hb=axis_ppm; spec_expt_hb=spec;

% Normalise the data
spec_expt_f= spec_expt_f/max(spec_expt_f);
spec_expt_ha=spec_expt_ha/max(spec_expt_ha);
spec_expt_hb=spec_expt_hb/max(spec_expt_hb);

% Set the guess
guess=[-14.4404    2.1802   10.0564   14.0478   38.0147   27.7945   18.4317  ...
         1.2295    3.1159  -15.1172   11.2721   49.6307   18.9158    4.5308  7.5739];

% Set optimiser options
options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf);

% Get the figure going
figure(); scale_figure([2.00 1.75]);

% Run the optimisation
answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f,  ...
                             axis_expt_ha,spec_expt_ha, ...
                             axis_expt_hb,spec_expt_hb,x),guess,options);

% Display the result
disp(answer);

end

% Least squares error function
function err=errfun(axis_expt_f, spec_expt_f,  ...
                    axis_expt_ha,spec_expt_ha, ...
                    axis_expt_hb,spec_expt_hb,params)

% Silence the output
sys.output='hush';
sys.disable={'hygiene'};

% Magnet induction 
sys.magnet=11.7464;

% Isotopes
sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C',     ...
              '1H',  '1H',  '19F', '1H',  '1H',  '1H',  'G',  'G', ...
              'G',   '1H',  '19F', '1H',  '1H',  'G',   'G',  'G'};  

% Chemical shifts
inter.zeeman.scalar=cell(1,23);
inter.zeeman.scalar{14}=    1.0092;
inter.zeeman.scalar{15}=    1.0092;
inter.zeeman.scalar{16}=    1.0092;
inter.zeeman.scalar{21}=    1.0092;
inter.zeeman.scalar{22}=    1.0092;
inter.zeeman.scalar{23}=    1.0092;
inter.zeeman.scalar{11}=    4.6834; 
inter.zeeman.scalar{17}=    4.6834;
inter.zeeman.scalar{10}= -184.1865;
inter.zeeman.scalar{18}= -184.1865;
inter.zeeman.scalar{8}=     1.7970;
inter.zeeman.scalar{9}=     1.7970;
inter.zeeman.scalar{13}=    1.6942;
inter.zeeman.scalar{20}=    1.6942;
inter.zeeman.scalar{19}=    1.6370;
inter.zeeman.scalar{12}=    1.6370;

% J-couplings
inter.coupling.scalar=cell(23);  
inter.coupling.scalar{19, 20}= params(1); 
inter.coupling.scalar{17, 8}=  params(2);
inter.coupling.scalar{17, 9}=  params(3);
inter.coupling.scalar{17, 19}= params(14);
inter.coupling.scalar{17, 20}= params(15);
inter.coupling.scalar{12, 13}= params(1);
inter.coupling.scalar{11, 8}=  params(3);
inter.coupling.scalar{11, 9}=  params(2);
inter.coupling.scalar{11, 12}= params(14);
inter.coupling.scalar{11, 13}= params(15);
inter.coupling.scalar{8,  9}=  params(10);
inter.coupling.scalar{10, 8}=  params(4); 
inter.coupling.scalar{10, 9}=  params(5); 
inter.coupling.scalar{10, 11}= params(12);
inter.coupling.scalar{10, 13}= params(7); 
inter.coupling.scalar{10, 12}= params(6);
inter.coupling.scalar{18, 9}=  params(4);
inter.coupling.scalar{18, 8}=  params(5);
inter.coupling.scalar{18, 17}= params(12);
inter.coupling.scalar{18, 20}= params(7); 
inter.coupling.scalar{18, 19}= params(6); 
inter.coupling.scalar{10, 18}= params(8);
inter.coupling.scalar{12,14}=  7.45;
inter.coupling.scalar{12,15}=  7.45;
inter.coupling.scalar{12,16}=  7.45;
inter.coupling.scalar{13,14}=  7.45;
inter.coupling.scalar{13,15}=  7.45;
inter.coupling.scalar{13,16}=  7.45;
inter.coupling.scalar{19,21}=  7.45;
inter.coupling.scalar{19,22}=  7.45;
inter.coupling.scalar{19,23}=  7.45;
inter.coupling.scalar{20,21}=  7.45;
inter.coupling.scalar{20,22}=  7.45;
inter.coupling.scalar{20,23}=  7.45;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas); 

% Sequence parameters (19F)
parameters_f.spins={'19F'};
parameters_f.rho0=state(spin_system,'L+','19F');
parameters_f.coil=state(spin_system,'L+','19F');
parameters_f.decouple={};
parameters_f.offset=-86700;
parameters_f.sweep=300;
parameters_f.npoints=512;
parameters_f.zerofill=2048;
parameters_f.axis_units='ppm';
parameters_f.invert_axis=1;

% Sequence parameters (1H A)
parameters_ha.spins={'1H'};
parameters_ha.rho0=state(spin_system,'L+',[11 17]);
parameters_ha.coil=state(spin_system,'L+',[11 17]);
parameters_ha.decouple={};
parameters_ha.offset=2335;
parameters_ha.sweep=128;
parameters_ha.npoints=256;
parameters_ha.zerofill=1024;
parameters_ha.axis_units='ppm';
parameters_ha.invert_axis=1;

% Sequence parameters (1H B)
parameters_hb.spins={'1H'};
parameters_hb.rho0=state(spin_system,'L+',[8 9]);
parameters_hb.coil=state(spin_system,'L+',[8 9]);
parameters_hb.decouple={};
parameters_hb.offset=890;
parameters_hb.sweep=128;
parameters_hb.npoints=256;
parameters_hb.zerofill=1024;
parameters_hb.axis_units='ppm';
parameters_hb.invert_axis=1;

% Simulation
fid_f= liquid(spin_system,@acquire,parameters_f, 'nmr');
fid_ha=liquid(spin_system,@acquire,parameters_ha,'nmr');
fid_hb=liquid(spin_system,@acquire,parameters_hb,'nmr');

% Apodisation and scaling
fid_f= params(9)* apodization(fid_f, 'exp-1d',6.0)/1e4;
fid_ha=params(13)*apodization(fid_ha,'gaussian-1d',10.0)/1e5;
fid_hb=params(11)*apodization(fid_hb,'gaussian-1d',6.5)/1e5;

% Fourier transform
spec_theo_f= real(fftshift(fft(fid_f,parameters_f.zerofill)));   spec_theo_f= spec_theo_f(end:-1:1);
spec_theo_ha=real(fftshift(fft(fid_ha,parameters_ha.zerofill))); spec_theo_ha=spec_theo_ha(end:-1:1);
spec_theo_hb=real(fftshift(fft(fid_hb,parameters_hb.zerofill))); spec_theo_hb=spec_theo_hb(end:-1:1);

% Axis unification and interpolation
axis_theo_f= sweep2ticks(parameters_f.offset, parameters_f.sweep, parameters_f.zerofill);
axis_theo_ha=sweep2ticks(parameters_ha.offset,parameters_ha.sweep,parameters_ha.zerofill);
axis_theo_hb=sweep2ticks(parameters_hb.offset,parameters_hb.sweep,parameters_hb.zerofill);
axis_theo_f= -2*pi*1e6*axis_theo_f/spin_system.inter.basefrqs(10);
axis_theo_ha=-2*pi*1e6*axis_theo_ha/spin_system.inter.basefrqs(9);
axis_theo_hb=-2*pi*1e6*axis_theo_hb/spin_system.inter.basefrqs(9);
spec_theo_f= interp1(axis_theo_f, spec_theo_f, axis_expt_f, 'pchip'); axis_theo_f= axis_expt_f;
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
err=  norm(spec_expt_f-spec_theo_f)^2+   ...
    2*norm(spec_expt_ha-spec_theo_ha)^2+ ...
    2*norm(spec_expt_hb-spec_theo_hb)^2;

end

