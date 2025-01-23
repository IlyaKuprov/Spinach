% Probability density sloshing around a harmonic oscillator
% in the presence of a gravitational pull twards the left.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function slosh_test_2()

% Set oscillator parameters
parameters.frc_cnst=2e3;  % N/m
parameters.par_mass=1.00; % kg
parameters.box_size=2.00; % m
parameters.n_points=100;  % points
parameters.grv_cnst=250;  % m/s^2

% Get the Hamiltonian
[H,~,xgrid]=oscillator(parameters);

% Get the initial state
psi=exp(-50*(xgrid-0.1).^2);

% Get the propagator
P=expm(full(-1i*H*0.001));

% Run the evolution
figure();
for n=1:1000
    plot(xgrid,abs(psi)+xgrid.^2,'r-'); hold on;
    plot(xgrid,xgrid.^2,'b-'); axis([-1 1 0 1.5]);
    kxlabel('particle coordinate, m');
    kylabel('probability density, a.u.'); kgrid;
    drawnow; hold off; pause(0.001); psi=P*psi;
end

end

