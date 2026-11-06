% Counterrotating field, Bloch-Siegert shift, and Floquet micromotion.
% Related to Section 4.1, Eq. (36), and Fig. 4 of Bock, Ernst,
% and Uhrig, https://doi.org/10.5194/mr-2026-10. This is a rotating-
% frame time-domain illustration, not a reproduction of their spectra.
% Angular frequencies are rad/s. The Fourier frequency is twice the
% resonant carrier frequency; order 2 means third order in the paper.
%
% Local driven_rhs callback: t is seconds; u is a 4-by-1 column-stacked
% propagator; H is a 1-by-3 signed-harmonic cell in rad/s; omega is
% rad/s; du is its 4-by-1 time derivative in inverse seconds.
%
% talos@spindynamics.org

function fme_bloch_siegert()

% One spin in Hilbert space
sys.isotopes={'1H'}; sys.magnet=0;
sys.parallel={'processes',1}; sys.scratch=tempdir;
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
Lx=full(operator(spin_system,'Lx',1));
Lz=full(operator(spin_system,'Lz',1));
Lp=full(operator(spin_system,'L+',1));

% Resonant rotating frame, including the counterrotating field
omega0=2*pi*10000; omega1=2*pi*2000; omega=2*omega0;
H={omega1*Lp'/2,omega1*Lx,omega1*Lp/2};
[F,K]=floquet_magnus(H,omega,2,0,'van_vleck');
fprintf('Leading Bloch-Siegert coefficient: %.6g Hz along Lz.\n',...
        omega1^2/(4*omega0)/(2*pi));

% Direct time ordering: integrate the two-by-two propagator
period=2*pi/omega; times=linspace(0,40*period,1601);
options=odeset('RelTol',1e-10,'AbsTol',1e-12);
[~,trajectory]=ode113(@(t,u)driven_rhs(t,u,H,omega),...
                      times,reshape(eye(2),[],1),options);

% Compare complete evolution with omission of both endpoint kicks
rho0=Lz; signal=zeros(numel(times),3);
K=cellfun(@full,K,'UniformOutput',false); kicks=cat(3,K{:});
harmonics=-(numel(K)-1)/2:(numel(K)-1)/2;
K0=sum(kicks,3);
for n=1:numel(times)
    Kt=sum(kicks.*reshape(exp(1i*harmonics*omega*times(n)),1,1,[]),3);
    Uref=reshape(trajectory(n,:),2,2);
    Ubare=expm(-1i*F*times(n));
    Ufull=expm(-1i*Kt)*Ubare*expm(1i*K0);
    signal(n,1)=real(trace(Lz*Uref*rho0*Uref'))/trace(Lz^2);
    signal(n,2)=real(trace(Lz*Ufull*rho0*Ufull'))/trace(Lz^2);
    signal(n,3)=real(trace(Lz*Ubare*rho0*Ubare'))/trace(Lz^2);
end
fprintf('Maximum signal error: complete %.3g; missing kicks %.3g.\n',...
        max(abs(signal(:,2)-signal(:,1))),...
        max(abs(signal(:,3)-signal(:,1))));

% Rotor-period analogues: endpoints and within-period micromotion
kfigure(); subplot(2,1,1);
plot(times/period,signal); kgrid;
kxlabel('Counterrotating-field periods'); kylabel('Normalised signal');
legend('Direct','Full FVV','No kicks','Location','southoutside',...
       'Orientation','horizontal','Interpreter','none','FontSize',9);
subplot(2,1,2);
plot(times/period,signal(:,2:3)-signal(:,1)); hold on;
plot(times(1:40:end)/period,signal(1:40:end,2)-signal(1:40:end,1),'ko');
kgrid; kxlabel('Counterrotating-field periods'); kylabel('Signal error');
legend('Full FVV','No kicks','Period boundaries','Location','southoutside',...
       'Orientation','horizontal');

end

% ODE callback: t [s], u [4x1], H {1x3} [rad/s], omega [rad/s]; du [4x1/s]
function du=driven_rhs(t,u,H,omega)

% Evaluate the Hamiltonian and Schrodinger equation
Ht=H{1}*exp(-1i*omega*t)+H{2}+H{3}*exp(1i*omega*t);
du=reshape(-1i*Ht*reshape(u,2,2),[],1);

end


