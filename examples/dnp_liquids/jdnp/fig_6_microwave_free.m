% A demonstration of Maria Grazia Concilio's microwave-free JDNP
% effect where a field ramp in combination with unequal relaxati-
% on rates of singlet-alpha and singlet-beta product states crea-
% tes nuclear magnetisation enhancement beyond the Boltzmann le-
% vel at both the starting and the final field.
%
% Calculation time: minutes
%
% maria-grazia.concilio@weizmann.ac.il
% i.kuprov@soton.ac.uk

function fig_6_microwave_free()

% Load the spin system
[sys,inter,bas]=system_specification();

% Magnet fields
start_field=14.09;
match_field=11.74;
final_field=9.39;

% Match exchange coupling to the midpoint field
inter.coupling.scalar{2,3}=match_field*(spin('E')+spin('1H'))/(2*pi);  
    
% Increase viscosity
inter.tau_c={2.2e-9};

% Get thermal equilibrium at starting field
sys.magnet=start_field;
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'labframe');
H0=hamiltonian(spin_system,'left');
rho_eq=equilibrium(spin_system,H0);

% Set up a field ramp and time step
field_grid=linspace(start_field,final_field,211); dt=1e-4;

% Preallocate the trajectory
traj=zeros(numel(rho_eq),numel(field_grid)+1,'like',1i);

% Set starting point
traj(:,1)=rho_eq;

% Run through the field ramp
for n=1:numel(field_grid)

    % Set the magneti field
    sys.magnet=field_grid(n);

    % Run the housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Get the generators
    H=hamiltonian(assume(spin_system,'labframe'));
    R=relaxation(spin_system);

    % Run the time evolution
    traj(:,n+1)=evolution(spin_system,H+1i*R,[],traj(:,n),...
                          dt,1,'final');          
 
end 

% Build component operators
unit=unit_state(spin_system);
Nz=state(spin_system,{'Lz'},{1});
E1z=state(spin_system,{'Lz'},{2});
E2z=state(spin_system,{'Lz'},{3});
NzE1z=state(spin_system,{'Lz','Lz'},{1,2});
NzE2z=state(spin_system,{'Lz','Lz'},{1,3});
LzE2z=state(spin_system,{'Lz','Lz'},{2,3});
E1mE2p=state(spin_system,{'L-','L+'},{2,3});
E1pE2m=state(spin_system,{'L+','L-'},{2,3});
NzE1zE2z=state(spin_system,{'Lz','Lz','Lz'},{1,2,3});
NzE1mE2p=state(spin_system,{'Lz','L-','L+'},{1,2,3});
NzE1pE2m=state(spin_system,{'Lz','L+','L-'},{1,2,3});

% Build alpha singlet and triplet states
Sa=(unit/8)+(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z;
Tpa=(unit/8)+(Nz/4)+(E1z/4)+(E2z/4)+(LzE2z/2)+(NzE2z/2)+(NzE1z/2)+NzE1zE2z;
T0a=(unit/8)+(Nz/4)+(E1mE2p/4)+(E1pE2m/4)-(LzE2z/2)+(NzE1mE2p/2)+(NzE1pE2m/2)-NzE1zE2z;
Tma=(unit/8)+(Nz/4)-(E1z/4)-(E2z/4)+(LzE2z/2)-(NzE2z/2)-(NzE1z/2)+NzE1zE2z;

% Build beta singlet and triplet states
Sb=(unit/8)-(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)+(NzE1mE2p/2)+(NzE1pE2m/2)+NzE1zE2z;
Tpb=(unit/8)-(Nz/4)+(E1z/4)+(E2z/4)+(LzE2z/2)-(NzE2z/2)-(NzE1z/2)-NzE1zE2z;
T0b=(unit/8)-(Nz/4)+(E1mE2p/4)+(E1pE2m/4)-(LzE2z/2)-(NzE1mE2p/2)-(NzE1pE2m/2)+NzE1zE2z;
Tmb=(unit/8)-(Nz/4)-(E1z/4)-(E2z/4)+(LzE2z/2)+(NzE2z/2)+(NzE1z/2)-NzE1zE2z;
    
% Build Nz-singlet and Nz-triplet states
SNz=(Nz/4)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z;
TpNz=(Nz/4)+(NzE1z/2)+(NzE2z/2)+NzE1zE2z;
T0Nz=(Nz/4)+(NzE1mE2p/2)+(NzE1pE2m/2)-NzE1zE2z;
TmNz=(Nz/4)-(NzE1z/2)-(NzE2z/2)+NzE1zE2z;

% Detection states
coils=[Tpa,Tpb,Tma,Tmb,T0a,T0b,Sa,Sb,SNz,TpNz,T0Nz,TmNz,E1z,E2z,Nz];

% Compute the observables
answer=real(coils'*traj);

% Plotting
figure(); subplot(1,3,1); scale_figure([1.75 0.5]);
time_axis=linspace(0,numel(field_grid)*dt,numel(field_grid)+1);
plot(time_axis,answer(1,:),'b-'); hold on;
plot(time_axis,answer(3,:),'b-');
plot(time_axis,answer(5,:),'b-');
plot(time_axis,answer(2,:),'r-');
plot(time_axis,answer(4,:),'r-');
plot(time_axis,answer(6,:),'r-'); kgrid;
xlabel('time / seconds','interpreter','latex');
ylabel('state population','interpreter','latex');
legend({'${T_{+,\alpha}}$','${T_{-,\alpha}}$','${T_{0,\alpha}}$',...
        '${T_{+,\beta}}$','${T_{-,\beta}}$','${T_{0,\beta}}$'},...
        'interpreter','latex','Location','northeast');
set(gca,'TickLabelInterpreter','latex'); axis tight;
subplot(1,3,2);
plot(time_axis,answer(7,:),'b-'); hold on;
plot(time_axis,answer(8,:),'r-'); kgrid;
xlabel('time / seconds','interpreter','latex');
ylabel('state population','interpreter','latex');
legend({'${S_{\alpha}}$','${S_{\beta}}$'},...
        'interpreter','latex','Location','southeast');
set(gca,'TickLabelInterpreter','latex'); axis tight;
subplot(1,3,3);
plot(time_axis,answer(15,:),'b-'); hold on;
xlabel('time / seconds','interpreter','latex');
ylabel('state population','interpreter','latex');
legend({'${N_{\rm{Z}}}$'},'interpreter','latex','Location','southeast');
set(gca,'TickLabelInterpreter','latex'); kgrid; axis tight;

end

