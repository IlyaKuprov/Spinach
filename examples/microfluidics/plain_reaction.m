% Non-linear reaction kinetics in a situation when there is
% no hydrodynamics, diffusion, or spin dynamics. This is in-
% tended as a stepping stone to the more complicated cases
% in the same directory of the Spinach example set.
%
% Calculation time: seconds.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function plain_reaction()

% No spin system here
spin_system=bootstrap();

% Rate constants, mol/(L*s)
k1=0.5;  % towards exo  
k2=0.1;  % towards endo

% Cycloaddition reaction generator, including solvent
K=@(t,x)(1i*[-k1*x(2)-k2*x(2)  0                0      0     0;      
              0               -k1*x(1)-k2*x(1)  0      0     0;           
              0                k1*x(1)          0      0     0;
              0                k2*x(1)          0      0     0; 
              0                0                0      0     0]);        
                     
% Kinetic time grid, 20 seconds
nsteps=200; tmax=20; dt=tmax/nsteps;
time_axis=linspace(0,tmax,nsteps+1); 
 
% Preallocate concentration trajectory
x=zeros(5,nsteps+1);

% Initial concentrations, mol/L
x(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]; 

% Concentration dynamics
for n=1:nsteps 
    x(:,n+1)=iserstep(spin_system,K,x(:,n),n*dt,dt,'LG4'); 
end

% Plot concentrations, excluding solvent
figure(); plot(time_axis,real(x(1:4,:))); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'cyclopentadiene','acrylonitrile', ...
         'endo-norbornene carbonitrile',...
         'exo-norbornene carbonitrile'}, ...
         'Location','northeast'); drawnow;

end

