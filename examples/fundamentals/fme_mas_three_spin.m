% 3-spin homonuclear dipolar MAS with complete AHT and Van Vleck evolution.
% Related to Sections 4.2-4.3, Eq. (39), and Figs. 9-11 of Bock, Ernst,
% and Uhrig, https://doi.org/10.5194/mr-2026-10. These time-domain
% examples use one crystallite and illustrative proton coordinates,
% not the paper's parameters or powder-averaged spectra. Order 2 below
% is third order in the paper. Both gauges include endpoint kicks
% so their intra-rotor evolution is compared on the same footing.
%
% States, propagation, and observables use native Spinach interfaces.
% The direct reference is checked by three-point step quadrature refinement.
%
% talos@spindynamics.org

function fme_mas_three_spin()

% Homonuclear spins, chemical shifts in ppm, coordinates in Angstrom
sys.isotopes={'1H','1H','1H'}; sys.magnet=9.4;
sys.parallel={'processes',1}; sys.scratch=tempdir; sys.output='hush';
inter.zeeman.scalar={0,2,5};
inter.coordinates={[0 0 0],[0 0 2.5],[2.2 0.6 1.0]};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Native rotor stack: rank-two interactions have five Fourier terms
parameters.axis=[sqrt(2/3) 0 1/sqrt(3)];
parameters.offset=0; parameters.spins={'1H'};
parameters.max_rank=2; parameters.rframes={};
parameters.orientation=[0.3 0.7 0.2]; parameters.masframe='magnet';
stack=rotor_stack(spin_system,parameters,'nmr');
stack=cellfun(@full,stack,'UniformOutput',false);
fourier_h=fftshift(fft(cat(3,stack{:}),[],3),3)/numel(stack);
H=cell(1,5);
for n=1:5
    H{n}=fourier_h(:,:,n);
end

% Use the native rotor sign and expand spin blocks, not the enlarged FP matrix
rate=25000; omega=-2*pi*rate; period=1/rate;
[Fa,Ka]=floquet_magnus(H,omega,2,0,'aht');
[Fv,Kv]=floquet_magnus(H,omega,2,0,'van_vleck');
rho0=state(spin_system,'Lx','all');
coil=coil_state(spin_system,'Lx','all'); normalisation=hdot(coil,rho0);
times=linspace(0,20*period,801); timestep=times(2)-times(1);
nsteps=numel(times)-1; input_harmonics=-2:2;

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

% Native effective trajectories retain the initial kicks in both gauges
Ka=cellfun(@full,Ka,'UniformOutput',false); kicks_a=cat(3,Ka{:});
Kv=cellfun(@full,Kv,'UniformOutput',false); kicks_v=cat(3,Kv{:});
harmonics=-(numel(Kv)-1)/2:(numel(Kv)-1)/2;
initial_a=step(spin_system,sum(kicks_a,3),rho0,-1);
initial_v=step(spin_system,sum(kicks_v,3),rho0,-1);
trajectory_a=evolution(spin_system,Fa,[],initial_a,timestep,nsteps,'trajectory');
trajectory_v=evolution(spin_system,Fv,[],initial_v,timestep,nsteps,'trajectory');
signal=zeros(numel(times),3); signal(:,1)=reference(:,2);

% Native endpoint transformations and Frobenius observable readout
for n=1:numel(times)
    phase=reshape(exp(1i*harmonics*omega*times(n)),1,1,[]);
    kick_a=sum(kicks_a.*phase,3); kick_v=sum(kicks_v.*phase,3);
    rho_a=step(spin_system,kick_a,trajectory_a{n},1);
    rho_v=step(spin_system,kick_v,trajectory_v{n},1);
    signal(n,2)=real(hdot(coil,rho_a))/normalisation;
    signal(n,3)=real(hdot(coil,rho_v))/normalisation;
end
fprintf('3-spin maximum signal error: AHT %.3g; complete FVV %.3g.\n',...
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


