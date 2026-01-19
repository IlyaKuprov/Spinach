% Nuclear spin Orientation via Electron spin Locking (NOVEL) and pulsed
% solid effect (SE), steady-state version. For futher information see:
%
%              https://doi.org/10.1016/0022-2364(88)90190-4
%                  https://doi.org/10.1063/1.5000528
% 
% Syntax (call from powder context)::
%
%           dnp=noveldnp_steady(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     H - Hamiltonian matrix, received from 
%         context function
%
%     R - relaxation superoperator, received 
%         from context function, must be ther-
%         malised to some finite temperature
%
%     K - kinetics superoperator, received 
%         from context function
%
%     parameters.irr_powers   - microwave amplitude (aka electron
%                             nutation frequency), Hz 
%
%     parameters.coil         - detection state(s)
%
%     parameters.contact_dur  - contact time, seconds
%
%     parameters.shot_spacing - delay between microwave irradiation periods
%
%     parameters.flippulse    - 0: Solid Effect (no flip pulse)
%                               1: NOVEL (90-degree flip pulse)
%
%     parameters.flipback     - 0: NOVEL without flipback pulse
%                               1: NOVEL with flipback pulse
%
%     parameters.addshift     - shift to center the field profile
%
%     parameters.el_offs      - microwave resonance offsets
%
% Output:
%
%     dnp                     - steady state observable on the de-
%                               tection state vector as a function
%                               of microwave resonance offset
%
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function dnp=noveldnp_steady(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Build electron operators
Ex=operator(spin_system,'Lx','E');
Ey=operator(spin_system,'Ly','E');
Ez=operator(spin_system,'Lz','E');

% Assemble the Liouvillian
L=H+1i*R+1i*K;

% Loop over offsets
dnp=zeros(size(parameters.el_offs),'like',1i);
for n=1:numel(parameters.el_offs)

    % Add rotating frame offset terms to the Liouvillian
    L_curr=L+2*pi*(parameters.el_offs(n)+parameters.addshift)*Ez;

    if parameters.flippulse % NOVEL
    
        % Add microwave irradiation along X during the 90-degree flip pulse
        L_pulse=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(0)+Ey*sind(0));
        
        % Add microwave irradiation along -Y during the contact pulse
        L_contact=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270));
    
        % Compute and clean up NOVEL propagator
        P=propagator(spin_system,L_contact,parameters.contact_dur)*...
          propagator(spin_system,L_pulse,parameters.pulse_dur);
        P=clean_up(spin_system,P,spin_system.tols.prop_chop);
    
        if parameters.flipback

            % Add microwave irradiation along -X during the 90-degree flipback pulse
            L_flipback=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(180)+Ey*sind(180));

            % Compute and clean up NOVEL w flipback propagator
            P=propagator(spin_system,L_flipback,parameters.pulse_dur)*P;
            P=clean_up(spin_system,P,spin_system.tols.prop_chop);

         end

    else % SE
 
        % Add microwave irradiation during the contact pulse along -Y
        L_contact=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270));
    
        % Compute and clean up SE propagator
        P=propagator(spin_system,L_contact,parameters.contact_dur);
        P=clean_up(spin_system,P,spin_system.tols.prop_chop);
    
    end

    % Shot spacing delay
    P=propagator(spin_system,L_curr,parameters.shot_spacing)*P;
    P=clean_up(spin_system,P,spin_system.tols.prop_chop);

    % Compute the steady state
    rho=steady(spin_system,P,[],'newton');
   
    % Get the observable at the steady state
    dnp(:,n)=gather(parameters.coil'*rho);

end

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'irr_powers')
    error('electron Rabi frequency must be specified in parameters.irr_powers variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'contact_dur')
    error('the contact time must be specified in parameters.contact_dur variable.');
end
if ~isfield(parameters,'shot_spacing')
    error('the delay between microwave irradiation periods must be specified in parameters.shot_spacing variable.');
end
if ((parameters.flippulse~=1)&&(parameters.flippulse~=0))
    error('parameters.flippulse can only be 0 or 1.');
end
if ((parameters.flipback~=1)&&(parameters.flipback~=0))
    error('parameters.flipback can only be 0 or 1.');
end
if ~isfield(parameters,'addshift')
    error('a shift to center the field profile must be specified in parameters.addshift variable.');
end
if ~isfield(parameters,'el_offs')
    error('the microwave resonance offsets must be specified in parameters.el_offs variable.');
end
end

% Nothing was your own except the few
% cubic centimetres inside your skull.
%
% George Orwell, "1984"

