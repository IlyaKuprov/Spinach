% Simultaneous fitting of multiple 51V MAS NMR spectra with 
% respect to the chemical shielding anisotropy and quadrupole
% coupling tensor parameters.
%
% Calculation time: hours, much faster with a Tesla A100 GPU.
%
% m.carravetta@soton.ac.uk
% i.kuprov@soton.ac.uk

function vfit_4sim_ik()

% Load and filter the data
s29=load('v12_29_dec15.spc'); s29=sgolayfilt(s29,3,51);  
s31=load('v12_31_dec15.spc'); s31=sgolayfilt(s31,3,51); 
s33=load('v12_33_dec15.spc'); s33=sgolayfilt(s33,3,51); 
s35=load('v12_35_dec15.spc'); s35=sgolayfilt(s35,3,51); 

% Set spectral ranges
kstart     = 6200;
kend_a     = 10000;
kend_b     = kstart+4096-1;
sp_range   = kstart:kend_a;
myrangeb   = kstart:kend_b;
nrpoints_a = kend_a-kstart+1;
nrpoints_b = 4096;

% Preprocess the spectra
S35=zeros(nrpoints_b,1); A35=s35(myrangeb,1);
S33=zeros(nrpoints_b,1); A33=s33(myrangeb,1);
S31=zeros(nrpoints_b,1); A31=s31(myrangeb,1);
S29=zeros(nrpoints_b,1); A29=s29(myrangeb,1);
S35(1:nrpoints_a,1)=s35(sp_range,2)/max(s35(sp_range,2)); 
S33(1:nrpoints_a,1)=s33(sp_range,2)/max(s33(sp_range,2)); 
S31(1:nrpoints_a,1)=s31(sp_range,2)/max(s31(sp_range,2)); 
S29(1:nrpoints_a,1)=s29(sp_range,2)/max(s29(sp_range,2)); 

% Set the initial guess
guess=[-669.0  564.0  0.255  82.0  180.0  19.0  3.72 0.62];

% Set optimiser options
options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf);

% Get a figure going
figure(); scale_figure([1.75 1.50]);

% Run the optimisation
fminsearch(@(x)errfun(A35,A33,A31,A29,S35,S33,S31,S29,x),guess,options);

end


% Least squares error function
function err=errfun(A35,A33,A31,A29,S35,S33,S31,S29,params)

% Display the parameters
disp(params);

% Silence Spinach
sys.output='hush';
sys.disable={'hygiene'};
% sys.enable={'gpu'};

% Data dimensions
sfO1 = 157.7815; nrpointsb = 4096;
swppm = A35(nrpointsb)-A35(1);
swhz = swppm * sfO1;

% Absorb parameters
CSiso=params(1);
CSaniso=-params(2);
CSeta=params(3);
CSAeul_a=params(4);
CSAeul_b=params(5);
CSAeul_c=params(6);
Qcc=params(7)*1.0e6;
Qeta=params(8);
lw=13000.0;

% Spin system
sys.isotopes={'51V'};

% Magnet field
sys.magnet=14.1;

% Chemical shifts
CSAa=CSiso - 0.5*CSaniso*(CSeta+1.0);
CSAb=CSiso + 0.5*CSaniso*(CSeta-1.0);
CSAc=CSiso + CSaniso;
inter.zeeman.eigs={[CSAa CSAb CSAc]+456.818}; 
inter.zeeman.euler={[CSAeul_a CSAeul_b CSAeul_c]};

% Quadrupolar couplings
inter.coupling.matrix{1,1}=eeqq2nqi(Qcc,Qeta,3.5,[0 0 0]);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Parameters
parameters.axis=[1 1 1];
parameters.max_rank=30;
parameters.grid='rep_2ang_200pts_oct';
parameters.sweep=swhz;
parameters.npoints=nrpointsb;
parameters.zerofill=nrpointsb;
parameters.offset=0.5*(A31(nrpointsb)+A31(1));
parameters.spins={'51V'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','51V');
parameters.coil=state(spin_system,'L+','51V');

% Simulation A
parameters.rate=41000;
fida=singlerot(spin_system,@acquire,parameters,'nmr');
fida=apodization(fida,'gaussian-1d',lw);
sim_speca=fftshift(fft(fida, parameters.zerofill));
sim_speca=sim_speca/max(sim_speca);

% Simulation B
parameters.rate=38500;
fidb=singlerot(spin_system,@acquire,parameters,'nmr');
fidb=apodization(fidb,'gaussian-1d',lw);
sim_specb=fftshift(fft(fidb,parameters.zerofill));
sim_specb=sim_specb/max(sim_specb);

% Simulation C
parameters.rate=36000;
fidc=singlerot(spin_system,@acquire,parameters,'nmr');
fidc=apodization(fidc,'gaussian-1d',lw);
sim_specc=fftshift(fft(fidc,parameters.zerofill));
sim_specc=sim_specc/max(sim_specc);

% Simulation D
parameters.rate=34000;
fidd=singlerot(spin_system,@acquire,parameters,'nmr');
fidd=apodization(fidd,'gaussian-1d',lw);
sim_specd=fftshift(fft(fidd,parameters.zerofill));
sim_specd=sim_specd/max(sim_specd);

% Plotting
subplot(2,2,1); plot(A35,real(S35),'ro','MarkerSize',1);
hold on; plot(A35,real(sim_speca)); axis tight;
kgrid; kxlabel('chemical shift, ppm'); hold off;
subplot(2,2,2); plot(A33,real(S33),'ro','MarkerSize',1);
hold on; plot(A33,real(sim_specb)); axis tight;
kgrid; kxlabel('chemical shift, ppm'); hold off;
subplot(2,2,3); plot(A31,real(S31),'ro','MarkerSize',1); 
hold on; plot(A31,real(sim_specc)); axis tight;
kgrid; kxlabel('chemical shift, ppm'); hold off;
subplot(2,2,4); plot(A29,real(S29),'ro','MarkerSize',1); 
hold on; plot(A29,real(sim_specd)); axis tight;
kgrid; kxlabel('chemical shift, ppm'); hold off; drawnow();

% Error functional
err=norm(real(S35)-real(sim_speca))^2+...
    norm(real(S33)-real(sim_specb))^2+...
    norm(real(S31)-real(sim_specc))^2+...
    norm(real(S29)-real(sim_specd))^2;

end

