% Radiation damping of an inverted proton: the Bloembergen-Pound
% tanh recovery and sech transverse envelope without relaxation.
% https://doi.org/10.1103/PhysRev.95.8
%
% Ported from Spynach, copyright (c) 2026 Alexej Jerschow, MIT licence:
% https://github.com/ajerschow/spynach/blob/main/LICENSE
%
% Calculation time: seconds
%
% Alexej Jerschow (original demonstrations)
% talos@spindynamics.org

function rd_closed_form()

% Magnet field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H'};

% Chemical shifts
inter.zeeman.scalar={0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotation generators and detection states
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
coil=[state(spin_system,'Lx','1H') ...
      state(spin_system,'Ly','1H') ...
      state(spin_system,'Lz','1H')];
rho_eq=coil(:,3); mz_eq=real(rho_eq'*rho_eq);

% Initial tilt and radiation damping rate
theta=170*pi/180; rate=40;
dt=5e-4; nsteps=600; time=(0:nsteps)*dt;
rho=step(spin_system,Ly,rho_eq,theta);
G=@(t,rho)rate*(real(coil(:,2)'*rho)*Lx-...
                        real(coil(:,1)'*rho)*Ly)/mz_eq;

% Fourth-order Lie-group product quadrature
magn=zeros(3,nsteps+1); magn(:,1)=real(coil'*rho)/mz_eq;
for n=1:nsteps
    rho=iserstep(spin_system,{G,time(n),'LG4'},rho,dt);
    magn(:,n+1)=real(coil'*rho)/mz_eq;
end

% Bloembergen-Pound closed form
t0=log(tan(theta/2))/rate;
exact=[sech(rate*(time-t0)); zeros(size(time)); ...
       tanh(rate*(time-t0))];
fprintf('Maximum closed-form deviation: %g\n',max(abs(magn-exact),[],'all'));

% Normalised magnetisation traces
kfigure();
subplot(1,2,1); plot(time,magn);
xlabel('time, seconds'); ylabel('magnetisation / equilibrium');
legend('m_x','m_y','m_z','Location','best'); kgrid;
subplot(1,2,2);
plot(time,magn(3,:),time,hypot(magn(1,:),magn(2,:)),...
     time,exact(3,:),'--',time,exact(1,:),'--');
xlabel('time, seconds'); ylabel('magnetisation / equilibrium');
legend('m_z','m_{perp}','tanh','sech','Location','best'); kgrid;

end

