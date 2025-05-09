% Beam DNP experiment from:
%
%             https://doi.org/10.1126/sciadv.abq0536
%
% Syntax (call from powder context):
%
%        contact_curve=beamdnp(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     H - Hamiltonian matrix, received from context function
%
%     R - relaxation superoperator, received from context function
%
%     K - kinetics superoperator, received from context function
%
%     parameters.irr_powers - microwave amplitude (aka electron
%                             nutation frequency), Hz 
%
%     parameters.rho0       - initial state
%
%     parameters.coil       - detection state
%
%     parameters.pulse_dur  - two pulse durations, seconds
%
%     parameters.nloops     - number of BEAM DNP blocks
%
% Output:
%
%     contact_curve         - time dependence of the coil state
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
%
% <https://spindynamics.org/wiki/index.php?title=beamdnp.m>

function contact_curve=beamdnp(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Build electron control operators
Ep=operator(spin_system,'L+','E'); 
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Microwave irradiation along Y during the 90-degree flip pulse
L_pulse=L+2*pi*parameters.irr_powers*Ey;
    
% Calculate the length of the 90-degree flip pulse
pulse_dur=1/(4*parameters.irr_powers);

% Apply the 90-degree flip pulse
rho=evolution(spin_system,L_pulse,[],parameters.rho0,pulse_dur,1,'final');

% Make repeated pulse Liouvillians
L1=L+2*pi*parameters.irr_powers*Ex;
L2=L-2*pi*parameters.irr_powers*Ex;

% Preallocate contact curve
contact_curve=zeros(1,parameters.nloops+1);
contact_curve(1)=hdot(parameters.coil,parameters.rho0);

% Adapt to the formalism
switch spin_system.bas.formalism

    % Hilbert space
    case 'zeeman-hilb'

        % Precompute the complete loop propagator
        P=propagator(spin_system,L2,parameters.pulse_dur(2))*...
          propagator(spin_system,L1,parameters.pulse_dur(1));

        % Run the sequence
        for k=1:parameters.nloops

            % Run a loop and record the observable
            rho=P*rho*P'; contact_curve(1+k)=hdot(parameters.coil,rho);

        end

    % Liouville space
    case {'zeeman-liouv','sphten-liouv'}

        % Precompute individual event propagators
        P1=propagator(spin_system,L1,parameters.pulse_dur(1));
        P2=propagator(spin_system,L2,parameters.pulse_dur(2));

        % Run the sequence
        for k=1:parameters.nloops

            % Run a loop and record the observable
            rho=P2*(P1*rho); contact_curve(1+k)=hdot(parameters.coil,rho);

        end

    otherwise

        % Complain and bomb out
        error('unknown formalism.');

end

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse duration should be specified in parameters.pulse_dur variable.');
end
if ~isfield(parameters,'nloops')
    error('the number of loops must be specified in parameters.nloops variable.');
end
if (~isnumeric(parameters.nloops))||(numel(parameters.nloops)~=1)||...
   (~isreal(parameters.nloops))||(parameters.nloops<1)||(mod(parameters.nloops,1)~=0)
    error('parameters.nloops should be a positive integer.');
end
end

% All we know about the new economic world tells 
% us that nations which train engineers will pre-
% vail over those which train lawyers. No nation
% has ever sued its way to greatness.
%
% Richard Lamm

