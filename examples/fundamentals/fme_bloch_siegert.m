% Counterrotating field, Bloch-Siegert shift, and Floquet micromotion.
% Related to Section 4.1, Eq. (36), and Fig. 4 of Bock, Ernst,
% and Uhrig, https://doi.org/10.5194/mr-2026-10. This is a rotating-
% frame time-domain illustration, not a reproduction of their spectra.
% Angular frequencies are rad/s. The Fourier frequency is twice the
% resonant carrier frequency; order 2 means third order in the paper.
%
% Time evolution uses evolution and step; hdot provides observables.
% The direct reference is checked by native three-point quadrature refinement.
%
% talos@spindynamics.org

function fme_bloch_siegert()

% One spin in Hilbert space
sys.isotopes={'1H'}; sys.magnet=0;
sys.parallel={'processes',1}; sys.scratch=tempdir; sys.output='hush';
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
Lx=full(operator(spin_system,'Lx',1));
Lp=full(operator(spin_system,'L+',1));

% Resonant rotating frame, including the counterrotating field
omega0=2*pi*10000; omega1=2*pi*2000; omega=2*omega0;
H={omega1*Lp'/2,omega1*Lx,omega1*Lp/2};
[F,K]=floquet_magnus(H,omega,2,0,'van_vleck');
fprintf('Leading Bloch-Siegert coefficient: %.6g Hz along Lz.\n',...
        omega1^2/(4*omega0)/(2*pi));

% Native states and detection operators keep the observable convention explicit
rho0=state(spin_system,'Lz',1); coil=coil_state(spin_system,'Lz',1);
normalisation=hdot(coil,rho0);
period=2*pi/omega; times=linspace(0,40*period,1601);
timestep=times(2)-times(1); nsteps=numel(times)-1; input_harmonics=-1:1;

% Assemble the physical Fourier series for native product quadrature
ham_data=cat(3,H{:});
ham_at=@(t)sum(ham_data.*reshape(exp(1i*input_harmonics*omega*t),1,1,[]),3);
reference=zeros(numel(times),2);
for pass=1:2

    % Refine the native fourth-order time ordering
    substeps=8*2^(pass-1); dt=timestep/substeps; rho=rho0;
    for n=1:numel(times)
        reference(n,pass)=real(hdot(coil,rho))/normalisation;
        if n==numel(times), break; end
        for k=1:substeps
            t=times(n)+(k-1)*dt;
            generators={ham_at(t),ham_at(t+dt/2),ham_at(t+dt)};
            rho=step(spin_system,generators,rho,dt);
        end
    end

end
ref_error=max(abs(reference(:,1)-reference(:,2)));
assert(ref_error<1e-8,'Native reference needs a finer time grid.');
fprintf('Native reference refinement difference: %.3g.\n',ref_error);

% Propagate the effective generator using the existing trajectory interface
K=cellfun(@full,K,'UniformOutput',false); kicks=cat(3,K{:});
harmonics=-(numel(K)-1)/2:(numel(K)-1)/2; K0=sum(kicks,3);
initial=step(spin_system,K0,rho0,-1);
trajectory=evolution(spin_system,F,[],initial,timestep,nsteps,'trajectory');
signal=zeros(numel(times),3); signal(:,1)=reference(:,2);
signal(:,3)=real(evolution(spin_system,F,coil,rho0,timestep,nsteps,...
                          'observable'))/normalisation;

% Apply endpoint kicks through step and read observables through hdot
for n=1:numel(times)
    Kt=sum(kicks.*reshape(exp(1i*harmonics*omega*times(n)),1,1,[]),3);
    rho=step(spin_system,Kt,trajectory{n},1);
    signal(n,2)=real(hdot(coil,rho))/normalisation;
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


