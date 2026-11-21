% Radiation damping of a weak heteronuclear doublet and a strongly
% coupled proton pair, compared with independent spectral lines
% sharing the same coil. The line picture is exact for the AX
% doublet, but does not reproduce the coupled AB dynamics.
%
% Ported from Spynach, copyright (c) 2026 Alexej Jerschow, MIT licence:
% https://github.com/ajerschow/spynach/blob/main/LICENSE
%
% Calculation time: seconds
%
% Alexej Jerschow (original demonstrations)
% talos@spindynamics.org

function rd_coupled_pair()

% Interaction and acquisition parameters
J=10; dnu=12; rate=40; theta=170*pi/180;
dt=5e-4; nsteps=600; time=(0:nsteps)*dt;
magn=zeros(3,nsteps+1,2,2);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Compare a weak AX doublet and a strongly coupled AB pair
for k=1:2

    % Magnet field
    sys.magnet=14.1;

    % Isotopes
    if k==1
        sys.isotopes={'1H','13C'};
    else
        sys.isotopes={'1H','1H'};
    end

    % Chemical shifts
    inter.zeeman.scalar={0 0};

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
    H=2*pi*J*operator(spin_system,{'Lz','Lz'},{1,2});
    if k==2
        H=H+pi*dnu*(operator(spin_system,'Lz',1)-...
                    operator(spin_system,'Lz',2))+...
            2*pi*J*(operator(spin_system,{'Lx','Lx'},{1,2})+...
                    operator(spin_system,{'Ly','Ly'},{1,2}));
    end

    % Quantum spin system with shared feedback
    G=@(t,rho)H+rate*(real(coil(:,2)'*rho)*Lx-...
                            real(coil(:,1)'*rho)*Ly)/mz_eq;
    rho=step(spin_system,Ly,rho_eq,theta);
    magn(:,1,k,1)=real(coil'*rho)/mz_eq;
    for n=1:nsteps
        rho=iserstep(spin_system,{G,time(n),'LG4'},rho,dt);
        magn(:,n+1,k,1)=real(coil'*rho)/mz_eq;
    end

    % Independent spectral lines and their equilibrium intensities
    if k==1
        offsets=[-J/2 J/2]; weights=[1 1]/2;
    else
        D=hypot(dnu,J);
        offsets=[-(D+J) -(D-J) D-J D+J]/2;
        weights=[1-J/D 1+J/D 1+J/D 1-J/D]/4;
    end
    % Single-spin system for the classical lines
    sys.isotopes={'1H'};
    inter.zeeman.scalar={0};
    line_system=create(sys,inter);
    line_system=basis(line_system,bas);
    Lx=operator(line_system,'Lx','1H');
    Ly=operator(line_system,'Ly','1H');
    Lz=operator(line_system,'Lz','1H');
    coil=[state(line_system,'Lx','1H') ...
          state(line_system,'Ly','1H') ...
          state(line_system,'Lz','1H')];
    rho_eq=coil(:,3); mz_eq=real(rho_eq'*rho_eq);
    niso=numel(offsets); detect=kron(weights,coil');
    Cx=kron(speye(niso),Lx); Cy=kron(speye(niso),Ly);
    H=kron(2*pi*spdiags(offsets',0,niso,niso),Lz);
    G=@(t,rho)H+rate*(real(detect(2,:)*rho)*Cx-...
                            real(detect(1,:)*rho)*Cy)/mz_eq;

    % Propagate the line ensemble through the same product quadrature
    rho=repmat(step(line_system,Ly,rho_eq,theta),niso,1);
    magn(:,1,k,2)=real(detect*rho)/mz_eq;
    for n=1:nsteps
        rho=iserstep(line_system,{G,time(n),'LG4'},rho,dt);
        magn(:,n+1,k,2)=real(detect*rho)/mz_eq;
    end
    fprintf('Case %d, maximum spin-system/line discrepancy: %g\n',...
            k,max(abs(magn(:,:,k,1)-magn(:,:,k,2)),[],'all'));
end

% Longitudinal recovery and transverse magnetisation
kfigure(); names={'AX doublet','AB pair'};
for k=1:2
    subplot(2,2,k);
    plot(time,magn(3,:,k,1),time,magn(3,:,k,2),'--');
    xlabel('time, seconds'); ylabel('m_z / equilibrium');
    title(names{k}); legend('spin system','lines','Location','best'); kgrid;
    subplot(2,2,k+2);
    plot(time,hypot(magn(1,:,k,1),magn(2,:,k,1)),...
         time,hypot(magn(1,:,k,2),magn(2,:,k,2)),'--');
    xlabel('time, seconds'); ylabel('m_{perp} / equilibrium');
    legend('spin system','lines','Location','northoutside'); kgrid;
end

end

