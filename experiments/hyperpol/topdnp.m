% Time-optimised pulsed DNP experiment from:
%
%                 https://doi.org/10.1126/sciadv.aav6909
% 
% Syntax (call from powder context):
%
%            contact_curve=topdnp(spin_system,parameters,H,R,K)
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
%     parameters.delay_dur  - delay_duration, seconds
%
%     parameters.nloops     - number of TOP DNP loops
%
% Output:
%
%     contact_curve         - time dependence of the coil state
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de
%
% <https://spindynamics.org/wiki/index.php?title=topdnp.m>

function contact_curve=topdnp(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Build electron control operators
Ep=operator(spin_system,'L+','E'); Ex=(Ep+Ep')/2; 

% Compose the Liouvillian
L=H+1i*R+1i*K;
 
% Add microwave irradiation along X
L_pulse=L+2*pi*parameters.irr_powers*Ex;

% Preallocate contact curve
contact_curve=zeros(1,parameters.nloops+1);
contact_curve(1)=hdot(parameters.coil,parameters.rho0);

% Adapt to the formalism
switch spin_system.bas.formalism

    % Hilbert space
    case 'zeeman-hilb'

        % Precompute TOP DNP loop propagator
        P=propagator(spin_system,L,parameters.delay_dur)*...
          propagator(spin_system,L_pulse,parameters.pulse_dur);

        % Run TOP DNP sequence
        rho=parameters.rho0;
        for k=1:parameters.nloops

            % Run a loop and record the observable
            rho=P*rho*P'; contact_curve(1+k)=hdot(parameters.coil,rho);

        end

    % Liouville space
    case {'zeeman-liouv','sphten-liouv'}

        % Run TOP DNP sequence
        rho=parameters.rho0;
        for k=1:parameters.nloops

            % Run a loop and record the observable
            rho=step(spin_system,L_pulse,rho,parameters.pulse_dur);
            rho=step(spin_system,L,rho,parameters.delay_dur);
            contact_curve(1+k)=hdot(parameters.coil,rho);

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
if ~isfield(parameters,'delay_dur')
    error('delay duration should be specified in parameters.pulse_dur variable.');
end
if ~isfield(parameters,'nloops')
    error('the number of loops must be specified in parameters.nloops variable.');
end
if (~isnumeric(parameters.nloops))||(numel(parameters.nloops)~=1)||...
   (~isreal(parameters.nloops))||(parameters.nloops<1)||(mod(parameters.nloops,1)~=0)
    error('parameters.nloops should be a positive integer.');
end
end

% Хотят ли русские войны?
% Спросите вы у тишины,
% Над ширью пашен и полей,
% И у берез, и тополей,
% Спросите вы у тех солдат,
% Что под березами лежат,
% И вам ответят их сыны,
% Хотят ли русские войны.
% 
% [...]
% 
% Да, мы умеем воевать,
% Но не хотим, чтобы опять
% Солдаты падали в бою
% На землю горькую свою.
% Спросите вы у матерей,
% Спросите у жены моей,
% И вы тогда понять должны
% Хотят ли русские войны.
% 
% Евгений Евтушенко, 1961

