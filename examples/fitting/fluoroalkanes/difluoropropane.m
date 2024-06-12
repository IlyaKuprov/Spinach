% Fitting of 1H and 19F NMR spectrums of 1,3-difluoropropane with 
% respect to J-couplings. See our paper for further details:
%
%           https://doi.org/doi/10.1021/acs.joc.4c00670
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% neil.wells@soton.ac.uk
% bruno.linclau@soton.ac.uk

function difluoropropane()

% Load experimental data
load('difluoropropane_fluorine.mat','axis_ppm','spec'); axis_expt_f=axis_ppm; spec_expt_f=spec;
load('difluoropropane_proton_a.mat','axis_ppm','spec'); axis_expt_ha=axis_ppm; spec_expt_ha=spec;
load('difluoropropane_proton_b.mat','axis_ppm','spec'); axis_expt_hb=axis_ppm; spec_expt_hb=spec;

% Normalise the data
spec_expt_f =spec_expt_f /max(spec_expt_f);
spec_expt_ha=spec_expt_ha/max(spec_expt_ha);
spec_expt_hb=spec_expt_hb/max(spec_expt_hb);

% Set the guess
guess=[0.3023    1.1605    0.7390   47.0061    5.7870   25.7705];

% Set optimiser options
options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf);

% Get the figure going
figure(); scale_figure([2.00 1.75]);

% Run the optimisation
answer=fminunc(@(x)errfun(axis_expt_f, spec_expt_f,...
                          axis_expt_ha,spec_expt_ha,...
                          axis_expt_hb,spec_expt_hb,x),guess,options);

% Display the result
disp(answer);

end

% Least squares error function
function err=errfun(axis_expt_f, spec_expt_f, ...
                    axis_expt_ha,spec_expt_ha,...
                    axis_expt_hb,spec_expt_hb,params)

% Silence the output
sys.output='hush';
sys.disable={'hygiene'};

% Magnet induction 
sys.magnet=11.7464;

% Isotopes
sys.isotopes={'1H','1H','19F','1H','1H','19F','1H','1H'};  

% Shifts
inter.zeeman.scalar=cell(1,8);
inter.zeeman.scalar{1}=4.6075;
inter.zeeman.scalar{2}=4.6075;
inter.zeeman.scalar{3}=-223.5314;
inter.zeeman.scalar{4}=2.1011;
inter.zeeman.scalar{5}=2.1011;
inter.zeeman.scalar{6}=-223.5314;
inter.zeeman.scalar{7}=4.6075;
inter.zeeman.scalar{8}=4.6075;

% J-couplings
inter.coupling.scalar=cell(8);  
inter.coupling.scalar{1,3}=params(4);
inter.coupling.scalar{2,3}=params(4);
inter.coupling.scalar{6,7}=params(4);
inter.coupling.scalar{6,8}=params(4);
inter.coupling.scalar{1,4}=params(5);
inter.coupling.scalar{1,5}=params(5);
inter.coupling.scalar{2,4}=params(5);
inter.coupling.scalar{2,5}=params(5);
inter.coupling.scalar{4,7}=params(5);
inter.coupling.scalar{4,8}=params(5);
inter.coupling.scalar{5,7}=params(5);
inter.coupling.scalar{5,8}=params(5);
inter.coupling.scalar{3,4}=params(6);
inter.coupling.scalar{3,5}=params(6);
inter.coupling.scalar{4,6}=params(6);
inter.coupling.scalar{5,6}=params(6);

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';
bas.sym_group={'S2','S2'};
bas.sym_spins={[1 2],[7 8]};

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas); 

% Sequence parameters (19F)
parameters_f.spins={'19F'};
parameters_f.rho0=state(spin_system,'L+',[3 6]);
parameters_f.coil=state(spin_system,'L+',[3 6]);
parameters_f.decouple={};
parameters_f.offset=-105059;
parameters_f.sweep=600;
parameters_f.npoints=2048;
parameters_f.zerofill=4096;
parameters_f.axis_units='ppm';
parameters_f.invert_axis=1;

% Sequence parameters (1H A)
parameters_ha.spins={'1H'};
parameters_ha.rho0=state(spin_system,'L+',[1 2 7 8]);
parameters_ha.coil=state(spin_system,'L+',[1 2 7 8]);
parameters_ha.decouple={};
parameters_ha.offset=2300;
parameters_ha.sweep=128;
parameters_ha.npoints=512;
parameters_ha.zerofill=1024;
parameters_ha.axis_units='ppm';
parameters_ha.invert_axis=1;

% Sequence parameters (1H B)
parameters_hb.spins={'1H'};
parameters_hb.rho0=state(spin_system,'L+',[4 5]);
parameters_hb.coil=state(spin_system,'L+',[4 5]);
parameters_hb.decouple={};
parameters_hb.offset=1050;
parameters_hb.sweep=300;
parameters_hb.npoints=1024;
parameters_hb.zerofill=2048;
parameters_hb.axis_units='ppm';
parameters_hb.invert_axis=1;

% Simulation
fid_f =liquid(spin_system,@acquire,parameters_f, 'nmr');
fid_ha=liquid(spin_system,@acquire,parameters_ha,'nmr');
fid_hb=liquid(spin_system,@acquire,parameters_hb,'nmr');

% Apodisation and scaling
fid_f= params(1)*apodization(fid_f, 'exp-1d',8.0)/3e2;
fid_ha=params(2)*apodization(fid_ha,'exp-1d',9.0)/3e2;
fid_hb=params(3)*apodization(fid_hb,'exp-1d',6.0)/3e2;

% Fourier transform
spec_theo_f= 0.05*real(fftshift(fft(fid_f, parameters_f.zerofill)));  spec_theo_f= spec_theo_f(end:-1:1);
spec_theo_ha=0.05*real(fftshift(fft(fid_ha,parameters_ha.zerofill))); spec_theo_ha=spec_theo_ha(end:-1:1);
spec_theo_hb=0.05*real(fftshift(fft(fid_hb,parameters_hb.zerofill))); spec_theo_hb=spec_theo_hb(end:-1:1);

% Axis unification and interpolation
axis_theo_f= sweep2ticks(parameters_f.offset, parameters_f.sweep, parameters_f.zerofill);
axis_theo_ha=sweep2ticks(parameters_ha.offset,parameters_ha.sweep,parameters_ha.zerofill);
axis_theo_hb=sweep2ticks(parameters_hb.offset,parameters_hb.sweep,parameters_hb.zerofill);
axis_theo_f= -2*pi*1e6*axis_theo_f /spin_system.inter.basefrqs(3);
axis_theo_ha=-2*pi*1e6*axis_theo_ha/spin_system.inter.basefrqs(1);
axis_theo_hb=-2*pi*1e6*axis_theo_hb/spin_system.inter.basefrqs(1);
spec_theo_f= interp1(axis_theo_f, spec_theo_f, axis_expt_f, 'pchip'); axis_theo_f= axis_expt_f;
spec_theo_ha=interp1(axis_theo_ha,spec_theo_ha,axis_expt_ha,'pchip'); axis_theo_ha=axis_expt_ha;
spec_theo_hb=interp1(axis_theo_hb,spec_theo_hb,axis_expt_hb,'pchip'); axis_theo_hb=axis_expt_hb;

% Plotting
subplot(3,1,1); plot(axis_expt_f,spec_expt_f,'r.'); hold on; plot(axis_theo_f,spec_theo_f,'b-'); 
hold off; xlim tight; set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm'); kgrid;
subplot(3,1,2); plot(axis_expt_ha,spec_expt_ha,'r.'); hold on; plot(axis_theo_ha,spec_theo_ha,'b-');
hold off; xlim tight; set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm'); kgrid;
subplot(3,1,3); plot(axis_expt_hb,spec_expt_hb,'r.'); hold on; plot(axis_theo_hb,spec_theo_hb,'b-');
hold off; xlim tight; set(gca,'XDir','reverse'); kxlabel('chemical shift, ppm'); kgrid;
disp(params); drawnow;

% Error functional
err=norm(spec_expt_f-spec_theo_f)^2+...
    norm(spec_expt_ha-spec_theo_ha)^2+...
    norm(spec_expt_hb-spec_theo_hb)^2;

end

