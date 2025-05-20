% TPPM DNP and its special case X-inverse-X (XiX) DNP experiment 
% from (https://doi.org/10.1021/jacs.1c09900). Syntax (call from
% powder context):
%
%        contact_curve=xixdnp(spin_system,parameters,H,R,K)
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
%     parameters.pulse_dur  - pulse duration, seconds
%
%     parameters.nloops     - number of XiX DNP blocks
%
%     parameters.phase      - phase of the second pulse
%
% Output:
%
%     contact_curve         - time dependence of the coil state
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de
%
% <https://spindynamics.org/wiki/index.php?title=xixdnp.m>

function contact_curve=xixdnp(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Build electron control operators
Ep=operator(spin_system,'L+','E'); 
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% Compose the Liouvillian
L=H+1i*R+1i*K;
    
% Add microwave irradiation along +X
L1=L+2*pi*parameters.irr_powers*Ex;

% Add microwave irradiation with user-specified phase
L2=L+2*pi*parameters.irr_powers*(Ex*cos(parameters.phase)+...
                                 Ey*sin(parameters.phase));

% Preallocate contact curve
contact_curve=zeros(1,parameters.nloops+1);
contact_curve(1)=hdot(parameters.coil,parameters.rho0);

% Grab initial condition and detection state
rho=parameters.rho0; coil=parameters.coil;

% Adapt to the formalism
switch spin_system.bas.formalism

    % Hilbert space
    case 'zeeman-hilb'

        % Precompute complete loop propagator
        P=propagator(spin_system,L2,parameters.pulse_dur)*...
          propagator(spin_system,L1,parameters.pulse_dur);

        % Run and record the observable
        for k=1:parameters.nloops
            rho=P*rho*P'; contact_curve(k+1)=hdot(coil,rho);
        end

    % Liouville space
    case {'zeeman-liouv','sphten-liouv'}

        % Precompute individual event propagators
        P1=propagator(spin_system,L1,parameters.pulse_dur);
        P2=propagator(spin_system,L2,parameters.pulse_dur);

        % Run and record the observable
        for k=1:parameters.nloops
            rho=P2*(P1*rho); contact_curve(k+1)=hdot(coil,rho);
        end

    otherwise

        % Complain and bomb out
        error('unknown formalism');

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

% Your PhD students will avenge me, Ilya.
%
% Peter Hore, to IK (who was then 
% a PhD student in his group)

