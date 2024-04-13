% A simple TOTAPOL based Cross Effect DNP system. Set to repro-
% duce Figure 2a from
%
%           http://dx.doi.org/10.1016/j.jmr.2011.09.047
%
% Intensity differences are due to a different relaxation model
% and minor inconsistencies between the stated geometry and the
% interaction amplitudes used in the original paper.
%
% Electron rotating frame simulation using Nottingham DNP rela-
% xation theory detailed in
%
%           http://dx.doi.org/10.1007/s00723-012-0367-0
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% alexander.karabanov@nottingham.ac.uk
% walter.kockenberger@nottingham.ac.uk

function cross_effect_freq_scan_2()

% Magnetic field
sys.magnet=3.4;

% Spin system
sys.isotopes={'E','E','1H'};
inter.zeeman.scalar={2.0023193 1.9981164 0.0000000};
inter.coordinates={[ 0.00   0.00   0.00];
                   [12.80   0.00   0.00];
                   [-3.12   0.00   3.12]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'nottingham'};
inter.rlx_keep='secular';
inter.equilibrium='zero';
inter.nott_r1e=1e2; 
inter.nott_r2e=1e5;
inter.nott_r1n=0.1; 
inter.nott_r2n=1e3;
inter.temperature=10;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.mw_pwr=2*pi*100e3;
parameters.mw_frq=2*pi*linspace(-350,350,5e4)*1e6;
parameters.coil=state(spin_system,'Lz','1H');
parameters.mw_oper=(operator(spin_system,'L-','E')+...
                    operator(spin_system,'L+','E'))/4;
parameters.ez_oper=operator(spin_system,'Lz','E');
parameters.orientation=[0 0 0];
parameters.method='lvn-backs';
parameters.needs={'aniso_eq'};
parameters.g_ref=inter.zeeman.scalar{1};

% Steady state simulation
answer=crystal(spin_system,@dnp_freq_scan,parameters,'esr');

% Plotting
figure(); plot(linspace(-350,350,5e4),real(answer)); kgrid;
axis tight; kxlabel('Microwave frequency offset, MHz');
kylabel('$S_\textrm{z}$ expectation value on $^{1}$H'); 
 
end

