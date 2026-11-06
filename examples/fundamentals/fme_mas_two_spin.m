% 2-spin homonuclear dipolar MAS with complete AHT and Van Vleck evolution.
% Related to Sections 4.2-4.3, Eq. (39), and Figs. 7-8 of Bock, Ernst,
% and Uhrig, https://doi.org/10.5194/mr-2026-10. These time-domain
% examples use one crystallite and illustrative proton coordinates,
% not the paper's parameters or powder-averaged spectra. Order 2 below
% is third order in the paper. Both gauges include endpoint kicks
% so their intra-rotor evolution is compared on the same footing.
%
% Local mas_rhs callback: t is seconds; u is a dim^2-by-1 column-stacked
% propagator; H is a 1-by-5 signed-harmonic cell in rad/s; omega is
% rad/s; dim is the Hilbert dimension; du is dim^2-by-1 in 1/s.
%
% talos@spindynamics.org

function fme_mas_two_spin()

% Homonuclear spins, chemical shifts in ppm, coordinates in Angstrom
sys.isotopes={'1H','1H'}; sys.magnet=9.4;
sys.parallel={'processes',1}; sys.scratch=tempdir;
inter.zeeman.scalar={0,2};
inter.coordinates={[0 0 0],[0 0 2.5]};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Native rotor stack: rank-two interactions have five Fourier terms
parameters.axis=[sqrt(2/3) 0 1/sqrt(3)];
parameters.offset=0; parameters.spins={'1H'};
parameters.max_rank=2; parameters.rframes={};
parameters.orientation=[0.3 0.7 0.2]; parameters.masframe='magnet';
[stack,phases]=rotor_stack(spin_system,parameters,'nmr');
H=cell(1,5);
for n=-2:2
    H{n+3}=zeros(size(stack{1}));
    for k=1:numel(phases)
        H{n+3}=H{n+3}+full(stack{k})*exp(-1i*n*phases(k))/numel(phases);
    end
end

% Use the native rotor sign and expand spin blocks, not the enlarged FP matrix
rate=25000; omega=-2*pi*rate; period=1/rate;
[Fa,Ka]=floquet_magnus(H,omega,2,0,'aht');
[Fv,Kv]=floquet_magnus(H,omega,2,0,'van_vleck');
rho0=full(operator(spin_system,'Lx','all')); coil=rho0;
dim=size(rho0,1); normalisation=trace(coil*rho0);

% Direct time ordering of the rank-two Hamiltonian, independent of F and K
times=linspace(0,20*period,801);
options=odeset('RelTol',1e-10,'AbsTol',1e-12);
[~,trajectory]=ode113(@(t,u)mas_rhs(t,u,H,omega,dim),...
                      times,reshape(eye(dim),[],1),options);

% Both intra-rotor points and rotor-synchronised points are retained
signal=zeros(numel(times),3);
Ka=cellfun(@full,Ka,'UniformOutput',false); kicks_a=cat(3,Ka{:});
Kv=cellfun(@full,Kv,'UniformOutput',false); kicks_v=cat(3,Kv{:});
harmonics=-(numel(Kv)-1)/2:(numel(Kv)-1)/2;
start_a=sum(kicks_a,3); start_v=sum(kicks_v,3);
for n=1:numel(times)
    Uref=reshape(trajectory(n,:),dim,dim);
    phase=reshape(exp(1i*harmonics*omega*times(n)),1,1,[]);
    kick_a=sum(kicks_a.*phase,3); kick_v=sum(kicks_v.*phase,3);
    Ua=expm(-1i*kick_a)*expm(-1i*Fa*times(n))*expm(1i*start_a);
    Uv=expm(-1i*kick_v)*expm(-1i*Fv*times(n))*expm(1i*start_v);
    signal(n,1)=real(trace(coil*Uref*rho0*Uref'))/normalisation;
    signal(n,2)=real(trace(coil*Ua*rho0*Ua'))/normalisation;
    signal(n,3)=real(trace(coil*Uv*rho0*Uv'))/normalisation;
end
fprintf('2-spin maximum signal error: AHT %.3g; complete FVV %.3g.\n',...
        max(abs(signal(:,2)-signal(:,1))),...
        max(abs(signal(:,3)-signal(:,1))));

% Compare complete intra-rotor evolution and mark stroboscopic samples
kfigure(); subplot(2,1,1);
plot(times/period,signal); kgrid;
kxlabel('Rotor periods'); kylabel('Normalised signal');
legend('Direct','Full AHT','Full FVV','Location','northeast',...
       'Interpreter','none','FontSize',9);
subplot(2,1,2);
plot(times/period,signal(:,2:3)-signal(:,1)); hold on;
plot(times(1:40:end)/period,signal(1:40:end,2:3)-signal(1:40:end,1),'o');
kgrid; kxlabel('Rotor periods'); kylabel('Signal error');
legend('AHT','FVV','AHT at nT','FVV at nT','Location','northeast',...
       'Interpreter','none','FontSize',9);

end

% ODE callback: t [s], u [dim^2x1], H {1x5} [rad/s], omega [rad/s]; du [1/s]
function du=mas_rhs(t,u,H,omega,dim)

% Sum the rank-two Fourier Hamiltonian and apply Schrodinger evolution
Ht=zeros(dim);
for n=-2:2
    Ht=Ht+H{n+3}*exp(1i*n*omega*t);
end
du=reshape(-1i*Ht*reshape(u,dim,dim),[],1);

end


