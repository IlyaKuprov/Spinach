% Hahn ESEEM of one spin-half nucleus coupled to a spin-one centre.
% The projected 0/+1 Hamiltonians give the analytic modulation in
% Eq. (29) of https://doi.org/10.1002/adts.202100254. Frequencies below
% are chosen model parameters, not a material-specific fit. Syntax:
%
%                     [answer,diagnostics]=cce_single_nucleus()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_single_nucleus()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E3','13C'};
sys.output='hush';
sys.parallel={'local',1};

% Nuclear Zeeman frequency and projected hyperfine components, rad/s
omega=-spin('13C')*sys.magnet;
apar=0.3*abs(omega);
aper=0.6*abs(omega);

% Standard Spinach tensors: electron g, nuclear ppm, and coupling Hz
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=[0 0 aper;0 0 0;0 0 apar]/(2*pi);

% Local Hilbert basis for retained clusters
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2;
parameters.transition.spins=1;
parameters.transition.ms=[0 1];
parameters.times=linspace(0,8*pi/abs(omega),101);
parameters.cce.order=1;
parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-9,'rtol',1e-9,...
                               'max_memory_bytes',2^28);

% Execute the named Hahn experiment through the cluster context
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Plot the measured complex echo amplitude
kfigure();
plot(1e6*parameters.times,real(answer),'b-',...
     1e6*parameters.times,imag(answer),'r--');
xlabel('total echo time / microseconds');
ylabel('normalised echo amplitude E(2τ)');
legend('real','imaginary');
kgrid();

end


