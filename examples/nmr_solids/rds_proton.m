% Proton-detected residual dipolar splitting of a 1H–14N pair under
% MAS, comparing the coupled system with its dipole-free reference.
% System and acquisition parameters: DOI 10.1016/j.jmro.2026.100218.
%
% talos@spindynamics.org

function rds_proton()

% Field and quadrupolar interaction
sys.magnet=2*pi*800e6/spin('1H');
sys.isotopes={'1H','14N'};
sys.parallel={'processes',64};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{2,2}=eeqq2nqi(3e6,0.8,1,[0 0 0]);

% Basis set and exact spin-space detection channel
bas.formalism='zeeman-liouv';
bas.approximation={'none'};
parameters.rate=50e3;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=15;
parameters.grid='rep_2ang_800pts_sph';
parameters.npoints=4096;
parameters.zerofill=8192;
parameters.sweep=20e3;
parameters.spins={'1H'};
parameters.offset=0;
parameters.axis_units='Hz';
parameters.rframes={{'1H',2},{'14N',2}};

% Compare the same acquisition with and without the signed dipolar coupling
dipole_hz=[-7000 0];
kfigure();
hold on
for n=1:numel(dipole_hz)

    % Spinach housekeeping
    inter.coupling.matrix{1,2}=dipole_hz(n)*diag([-1 -1 2]);
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Acquire proton signal with numerical rotating frames retaining mixed terms
    parameters.rho0=state(spin_system,'Lx',1);
    parameters.coil=coil_state(spin_system,'L+',1);
    fid=singlerot(spin_system,@acquire,parameters,'labframe');

    % Apply identical processing to both signals
    fid=apodisation(spin_system,fid,{{'exp',15}});
    spectrum=fftshift(fft(fid,parameters.zerofill));
    plot_1d(spin_system,real(spectrum),parameters);

end
klegend({'Dipole on','Dipole off'},'Location','Best');

end


