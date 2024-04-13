% A plot of a typical random walk on a sphere.
%
% jpresteg@uga.edu
% i.kuprov@soton.ac.uk

function rwalk_plot()

% tau_c and number of steps per tau_c
tau_c=6e-9; tc_steps=500; 
dt=tau_c/tc_steps;

% Total number of steps
tot_nsteps=5000; 

% Euler angles of random walk
rng('shuffle');
eulers=rwalk(tot_nsteps,tau_c,dt); 

% Trajectory preallocation
traj=zeros(3,size(eulers,1));

% Trajectory
for n=1:size(eulers,1)
    traj(:,n)=euler2dcm(eulers(n,:))*[0; 0; 1];
end

% Plotting
plot3(traj(1,:),traj(2,:),traj(3,:));
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
axis square; kgrid; box on;

end

