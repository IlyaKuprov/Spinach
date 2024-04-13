% Nuclear spin Orientation via Electron spin Locking (NOVEL) and pulsed
% solid effect (SE). For futher information see:
%
%              https://doi.org/10.1016/0022-2364(88)90190-4
%                  https://doi.org/10.1063/1.5000528
% 
% Syntax (call from powder context):
%
%          contact_curve=novel_se(spin_system,parameters,H,R,K)
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
%     parameters.timestep   - time step of the contact curve, s
%
%     parameters.nsteps     - number of time steps in the con-
%                             tact curve
%
%     parameters.flippulse  - 0: Solid Effect (no flip pulse)
%                             1: NOVEL (90-degree flip pulse) 
%
% Output:
%
%     contact_curve         - time dependence of the coil state
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
%
% <https://spindynamics.org/wiki/index.php?title=novel_se.m>

function contact_curve=novel_se(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Build electron control operators
Ep=operator(spin_system,'L+','E');
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% Compose the Liouvillian
L=H+1i*R+1i*K;

% 90-degree flip pulse
if parameters.flippulse
    
    % Add microwave irradiation along X during the 90-degree flip pulse
    L_pulse=L+2*pi*parameters.irr_powers*(Ex*cosd(0)+Ey*sind(0));
    
    % Calculate the length of the 90-degree flip pulse
    pulse_dur=1/(4*parameters.irr_powers);

    % Apply the 90-degree flip pulse
    rho=evolution(spin_system,L_pulse,[],parameters.rho0,pulse_dur,1,'final');

else
    
    % Do nothing
    rho=parameters.rho0;

end

% Set the spin lock direction to -Y
L_slock=L+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270));

% Compute the contact curve
contact_curve=evolution(spin_system,L_slock,parameters.coil,rho,...
                        parameters.timestep,parameters.nsteps,'observable');

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
if ~isfield(parameters,'nsteps')
    error('number of steps should be specified in parameters.nsteps variable.');
end
if (~isnumeric(parameters.nsteps))||(numel(parameters.nsteps)~=1)||...
   (~isreal(parameters.nsteps))||(parameters.nsteps<1)||...
   (mod(parameters.nsteps,1)~=0)
    error('parameters.nsteps should be a positive integer.');
end
if ((parameters.flippulse~=1)&&(parameters.flippulse~=0))
    error('parameters.flippulse can only be 0 or 1.');
end
end

% 10 Whatsoever thy hand findeth to do, do it with thy might; for there is
%    no work, nor device, nor knowledge, nor wisdom, in the grave, whither
%    thou goest.
%
% 11 I returned, and saw under the sun, that the race is not to the swift,
%    nor the battle to the strong, neither yet bread to the wise, nor yet
%    riches to men of understanding, nor yet favour to men of skill; but 
%    time and chance happeneth to them all.
%
% 12 For man also knoweth not his time: as the fishes that are taken in an
%    evil net, and as the birds that are caught in the snare; so are the 
%    sons of men snared in an evil time, when it falleth suddenly upon them.
%
% Ecclesiastes 9

