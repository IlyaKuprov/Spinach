% Projected Hahn echo of the two-nucleus silicon donor regression model.
% Syntax:
%
%             [answer,diagnostics]=cce_sip_small()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_sip_small()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','29Si','29Si'};
sys.output='hush';
sys.parallel={'local',1};
sys.tols.cce_atol=1e-8;
sys.tols.cce_rtol=1e-8;
sys.tols.cce_mem=2^28;

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[8 4 0];[12 2 4]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
radius=25;
valley_wavevector=0.85*2*pi/5.43;
for n=2:3

    % Hydrogenic envelope and six-valley interference, illustrative scale
    position=inter.coordinates{n};
    contact=2e5*exp(-2*norm(position)/radius)*...
            (sum(cos(valley_wavevector*position))/3)^2;

    % Retain observed Sz only for the explicitly projected donor model
    inter.coupling.matrix{1,n}=diag([0 0 contact]);

end

% Explicit dipolar coupling between the bath nuclei in Hz
inter.coupling.matrix{2,3}=diag([200 200 -400]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:3;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=2;
parameters.cce.approximation='projected';

% Execute the named Hahn echo through the cluster-expansion context
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end


