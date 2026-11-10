% Radiation damping of a Gaussian inhomogeneous proton line.
% All isochromats share one coil; their weighted transverse
% magnetisation generates the feedback field on every spin.
%
% Ported from Spynach (A. Jerschow, MIT licence).
% https://github.com/ajerschow/spynach
%
% Calculation time: seconds
%
% Alexej Jerschow (original demonstrations)
% Talos (MATLAB port)
% Copyright (c) 2026 Alexej Jerschow; MIT licence, see LICENSE.md.

function rd_inhomogeneous()

% Spin system and basis set
sys.magnet=14.1; sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv'; bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotation generators and detection states
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
Lz=operator(spin_system,'Lz','1H');
coil=[state(spin_system,'Lx','1H') ...
      state(spin_system,'Ly','1H') ...
      state(spin_system,'Lz','1H')];
rho_eq=coil(:,3); mz_eq=real(rho_eq'*rho_eq);

% Offset grid, linewidths, and acquisition parameters
offsets=linspace(-120,120,81); widths=[0 10 30 80];
rate=40; theta=170*pi/180; dt=5e-4; nsteps=600;
time=(0:nsteps)*dt; niso=numel(offsets);
rho0=step(spin_system,Ly,rho_eq,theta);
Cx=kron(speye(niso),Lx); Cy=kron(speye(niso),Ly);
H=kron(2*pi*spdiags(offsets',0,niso,niso),Lz);
magn=zeros(3,nsteps+1,numel(widths));

% Gaussian ensembles with the same equilibrium magnetisation
for k=1:numel(widths)
    if widths(k)==0
        weights=double(offsets==0);
    else
        weights=exp(-4*log(2)*(offsets/widths(k)).^2);
    end
    weights=weights/sum(weights);
    detect=kron(weights,coil');
    G=@(t,rho)H+rate*(real(detect(2,:)*rho)*Cx-...
                            real(detect(1,:)*rho)*Cy)/mz_eq;
    rho=repmat(rho0,niso,1);
    magn(:,1,k)=real(detect*rho)/mz_eq;
    for n=1:nsteps
        rho=iserstep(spin_system,{G,time(n),'LG4'},rho,dt);
        magn(:,n+1,k)=real(detect*rho)/mz_eq;
    end
end

% Dephasing weakens the shared radiation damping field
kfigure();
subplot(1,2,1); plot(time,squeeze(magn(3,:,:)));
xlabel('time, seconds'); ylabel('m_z / equilibrium');
legend('0 Hz','10 Hz','30 Hz','80 Hz','Location','best'); kgrid;
subplot(1,2,2); plot(time,squeeze(hypot(magn(1,:,:),magn(2,:,:))));
xlabel('time, seconds'); ylabel('m_{perp} / equilibrium');
legend('0 Hz','10 Hz','30 Hz','80 Hz','Location','best'); kgrid;

end

% The purpose of computing is insight, not numbers.
% Richard Hamming


