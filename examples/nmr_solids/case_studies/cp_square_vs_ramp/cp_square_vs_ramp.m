% 1H-15N cross-polarisation experiment in the doubly rotating
% frame using (a) fixed amplitude CP; (b) linearly ramped CP;
% (c) tangent-ramped CP. Static powder simulation demonstra-
% ting the advantages of ramped cross-polarisation. For fur-
% ther information, see:
%
%         https://doi.org/10.1016/0009-2614(94)00470-6
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function cp_square_vs_ramp()

% System specification
sys.magnet=9.394;
sys.isotopes={'15N','1H'};
          
% Interactions
inter.zeeman.scalar={0.00 0.00};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.05]};
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Common experiment parameters
parameters.time_steps=2e-6*ones(1,500);
parameters.irr_opers={operator(spin_system,'Ly','1H') ...
                      operator(spin_system,'Lx','15N')};
parameters.exc_opers={operator(spin_system,'Lx','1H')};
parameters.coil=state(spin_system,'Lx','15N');
parameters.grid='rep_2ang_6400pts_sph';
parameters.needs={'aniso_eq'};
parameters.spins={'15N'};

% Simulate fixed amplitude CP
irr_powers_a=[5e4*ones(1,500); 5e4*ones(1,500)];
parameters.irr_powers=irr_powers_a;
fid_a=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Simulate linearly ramped amplitude CP
ramp_up=linspace(0,1,500); ramp_down=fliplr(ramp_up);
irr_powers_b=[5e4*ramp_down; 5e4*ramp_up];
parameters.irr_powers=irr_powers_b;
fid_b=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Simulate tangent ramped amplitude CP
ramp_up=tan(linspace(-1.4,1.4,500));
ramp_up=ramp_up-min(ramp_up); 
ramp_up=ramp_up/max(ramp_up);
ramp_down=fliplr(ramp_up);
irr_powers_c=[5e4*ramp_down; 5e4*ramp_up];
parameters.irr_powers=irr_powers_c;
fid_c=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Plotting
figure(); scale_figure([1.5 0.75]);
time_axis=[0 cumsum(parameters.time_steps)];
subplot(1,2,1); plot(time_axis(2:end),[irr_powers_a;
                                       irr_powers_b;
                                       irr_powers_c]);
xlim tight; kgrid; kxlabel('time, seconds'); 
kylabel('nutation frequency, Hz'); ylim([-5e4 6e4]);
klegend({'constant ampl.','constant ampl.',...
         'linear ramp',   'linear ramp',   ...
         'tangent ramp',  'tangent ramp'},'Location','South');
subplot(1,2,2); plot(time_axis,real([fid_a; fid_b; fid_c]));
xlim tight; kgrid; kxlabel('time, seconds'); 
kylabel('$S_{\rm{X}}$ expectation value on $^{15}N$'); 
klegend({'constant ampl.','linear ramp','tangent ramp'},...
         'Location','SouthWest');

end

