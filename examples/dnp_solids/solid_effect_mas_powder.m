% A MAS DNP simulation performed as described in Fred Mentink-
% Vigier's paper (Spinach rotation conventions are different):
%
%         http://dx.doi.org/10.1016/j.jmr.2015.07.001
%
% Steady state DNP simulation for a powder.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function solid_effect_mas_powder()

% Magnet field
sys.magnet=9.403;

% Spin specification
sys.isotopes={'E','1H'};

% Interactions
inter.zeeman.eigs{1}=[2.00614  2.00194 2.00988];
inter.zeeman.euler{1}=pi*[253.6 105.1 123.8]/180;
inter.zeeman.eigs{2}=[0.00 0.00 0.00];
inter.zeeman.euler{2}=[0.00 0.00 0.00];
inter.coordinates={[0.00 0.00 0.00],...
                   [0.00 0.00 3.00]};

% Relaxation parameters
inter.relaxation={'weizmann'};
inter.weiz_r1e=1/0.3e-3;
inter.weiz_r1n=1/4.0;
inter.weiz_r2e=1/1.0e-6;
inter.weiz_r2n=1/0.2e-3;
inter.weiz_r1d=zeros(2,2);
inter.weiz_r2d=zeros(2,2);
inter.temperature=100;
inter.equilibrium='dibari';
inter.rlx_keep='secular';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Stack generation parameters
parameters.spins={'E'};
parameters.rate=12.5e3;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=800;
parameters.mw_pwr=2*pi*0.85e6;
parameters.mw_frq=-263.366e9;
parameters.mw_time=1.0;
parameters.grid='rep_2ang_100pts_sph';
parameters.coil=state(spin_system,'Lz','1H');
parameters.verbose=0;

% Run the MAS DNP simulation
answer=masdnp(spin_system,parameters);
disp(['Steady state DNP enhancement: ' num2str(answer)]);

end

