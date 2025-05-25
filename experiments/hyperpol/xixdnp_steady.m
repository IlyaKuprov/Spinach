% TPPM DNP and its special case X-inverse-X (XiX) DNP experiment 
% from (https://doi.org/10.1021/jacs.1c09900), steady state ver-
% sion. Syntax (call from powder context):
%
%        dnp=xixdnp_steady(spin_system,parameters,H,R,K)
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
%                               nutation frequency), Hz 
%
%     parameters.coil         - detection state vector
%
%     parameters.pulse_dur    - pulse duration, seconds
%
%     parameters.phase        - phase of each second pulse in radians
%
%     parameters.nloops       - number of XiX/TPPM DNP blocks, the
%                               calculation is faster when this is
%                               an integer power of 2
%
%     parameters.shot_spacing - delay between microwave irradiation
%                               periods, seconds
%
%     parameters.addshift     - shift of the centre of the field
%                               profile, Hz
%
%     parameters.el_offs      - microwave resonance offsets, a vector
%                               of frequencies in Hz
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
%
% <https://spindynamics.org/wiki/index.php?title=xixdnp_steady.m>

function dnp=xixdnp_steady(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Build electron control operators
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

    % Add microwave irradiation along +X
    L1=L_curr+2*pi*parameters.irr_powers*Ex;

    % Add microwave irradiation with user-specified phase
    L2=L_curr+2*pi*parameters.irr_powers*(Ex*cos(parameters.phase)+...
                                          Ey*sin(parameters.phase));

    % Compute the contact time block propagator
    P=propagator(spin_system,L2,parameters.pulse_dur)*...
      propagator(spin_system,L1,parameters.pulse_dur);
    P=clean_up(spin_system,P,spin_system.tols.prop_chop);

    % Send the problem to GPU if necessary
    if ismember('gpu',spin_system.sys.enable), P=gpuArray(P); end

    % Use propagator squaring to do 2^N blocks
    power_of_two=log2(parameters.nloops);
    if mod(power_of_two,1)==0
        for m=1:power_of_two
            P=clean_up(spin_system,P*P,spin_system.tols.prop_chop); 
        end
    else
        P=clean_up(spin_system,mpower(P,parameters.nloops),...
                               spin_system.tols.prop_chop);
    end

    % Get the problem from GPU if necessary
    if isa(P,'gpuArray'), P=gather(P); end

    % Add the delay after the contact time
    P=propagator(spin_system,L_curr,parameters.shot_spacing)*P;
    P=clean_up(spin_system,P,spin_system.tols.prop_chop);

    % Compute the stroboscopic steady state
    rho=steady(spin_system,P,[],[],'newton');
   
    % Get the observable at the steady state
    dnp(n)=gather(parameters.coil'*rho);

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
if ~isfield(parameters,'pulse_dur')
    error('pulse duration must be specified in parameters.pulse_dur variable.');
end
if ~isfield(parameters,'phase')
    error('second pulse phase must be specified in parameters.phase variable.');
end
if ~isfield(parameters,'nloops')
    error('nubmer of XiX/TPPM blocks must be specified in parameters.nloops variable.');
end
if ~isfield(parameters,'shot_spacing')
    error('delay between MW irradiation periods must be specified in parameters.shot_spacing variable.');
end
if parameters.shot_spacing<0
    error('parameters.shot_spacing must be a non-negative real number.');
end
if ~isfield(parameters,'addshift')
    error('field profile centre shift must be specified in parameters.addshift variable.');
end
if ~isfield(parameters,'el_offs')
    error('microwave resonance offsets must be specified in parameters.el_offs variable.');
end
end

% My foreign policy is to be able to take a 
% ticket at Victoria Station and go anywhere 
% I damn well please.
%
% Ernest Bevin, UK Foreign Secretary 1945-1951

