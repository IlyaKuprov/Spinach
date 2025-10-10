% Time dependence of LzSz and Lz+Sz spin orders in a para-
% hydrogen molecule coordinated to a nickel cage that cre-
% ates large chemical shift anisotropy. Done for Gloggler
% group, paper link coming in due course; see also the Ma-
% thematica worksheet.
%
% i.kuprov@soton.ac.uk

function rlx_trajectory()

% Magnet field, Tesla
sys.magnet=18.7893;

% Isotopes
sys.isotopes={'1H','1H'};

% Coordinates (for dipole tensor)
inter.coordinates={[ 0.1399 16.6491 21.2341];
                   [-1.1060 18.4823 21.4751]};

% Zeeman interaction tensors, traces subtracted
inter.zeeman.matrix={[61.897647  13.533565   1.607280
                      11.111659  47.381098 -25.571782
                      -2.442233 -20.906587  28.555729], ...
                     [39.335394  18.002885 -24.966389
                      20.627469  55.372869   5.145218
                     -20.402191   0.581529  43.347552]};
inter.zeeman.matrix=cellfun(@remtrace,inter.zeeman.matrix,...
                            'UniformOutput',false);

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={500e-12};

% Formalism and basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Assumptions
spin_system=assume(spin_system,'nmr');

% Build the Liouvillian
L=hamiltonian(spin_system)+1i*relaxation(spin_system);

% Initial state
rho=singlet(spin_system,1,2);

% Detection states
A=state(spin_system,{'L+','L-'},{1 2})+...
  state(spin_system,{'L-','L+'},{1 2});
B=state(spin_system,{'Lz','Lz'},{1 2});
C=state(spin_system,'Lz','all');

% Time evolution
traj=evolution(spin_system,L,[A/2 4*B C],rho,2e-3,1000,'multichannel');

% Plotting
time_axis=linspace(0,2,1001);
figure(); plot(time_axis,real(traj));
kxlabel('time, seconds'); kgrid;
kylabel('observable'); xlim tight; 
klegend({'$\bf{L}_{+}\bf{S}_{-}+\bf{L}_{-}\bf{S}_{+}$',...
         '$\bf{L}_\textrm{Z}\bf{S}_\textrm{Z}$',...
         '$\bf{L}_\textrm{Z}+\bf{S}_\textrm{Z}$'},...
         'Location','SouthEast','FontSize',14);
scale_figure([1.0 0.75]);
                   
end

