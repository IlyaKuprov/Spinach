% Davies ENDOR sequence with explicit soft pulses and all of the atten-
% dant effects, such as orientation selection. Soft pulses are simula-
% ted using the Fokker-Planck formalism. Syntax:
%
%           answer=endor_davies(spin_system,parameters,H,R,K)
%
% Parameters:
%
% The following parameters refer to the electron pi pulse. The duration
% of the electron pi/2 pulse is obtained by halving parameters.e_dur:
%
%     parameters.e_frq - frequency of the electron pulse, Hz
%
%     parameters.e_phi - phase of the electron pulse, rad
%
%     parameters.e_pwr - power of the electron pulse, rad/s
%
%     parameters.e_dur - duration of the electron pulse, s
%
%     parameters.e_rnk - Fokker-Planck cut-off rank for
%                        the electron pulse
%
% The following parameters refer to the nuclei pulse:
%
%     parameters.n_frq  - vector of frequencies for the nuclei
%                         pulse, in Hz. The answer is returned
%                         as a vector of the same dimension.
%
%     parameters.n_phi  - phase of the nuclei pulse, rad
%
%     parameters.n_pwr  - power of the nuclei pulse, rad/s
%
%     parameters.n_dur  - duration of the nuclei pulse, s
%
%     parameters.n_rnk  - Fokker-Planck cut-off rank for
%                         the nuclei pulse
%
%     parameters.method - method to use during the call
%                         to shaped_pulse_af()
%
%     H    - Hamiltonian matrix, received from context function
%
%     R    - relaxation superoperator, received from context function
%
%     K    - kinetics superoperator, received from context function
%
% Outputs:
%
%     answer - amplitude detected on the coil state for each 
%              frequency of the nuclear pulse
%
% Note: Fokker-Planck ranks should be increased until convergence is 
%       achieved in the output. The same applies to the size of the 
%       spherical grid.
%
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=endor_davies.m>

function answer=endor_davies(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Pulse operators
Ep=operator(spin_system,'L+',parameters.spins{1});
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% Pulse operators
Np=operator(spin_system,'L+',parameters.spins{2});
Nx=(Np+Np')/2; Ny=(Np-Np')/2i;

% Soft pulse frequencies
parameters.e_frq=parameters.e_frq-parameters.offset(1);
parameters.n_frq=parameters.n_frq+spin(parameters.spins{2})*spin_system.inter.magnet/(2*pi);

% Pi pulse on the electron
rho0=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0, parameters.e_frq,...
                                         parameters.e_pwr,parameters.e_dur,...
                                         parameters.e_phi,parameters.e_rnk,...
                                         parameters.method);

% Preallocate the answer
answer=zeros(size(parameters.n_frq),'like',1i);

% Loop over radiofrequency offsets
parfor n=1:numel(parameters.n_frq)
             
    % Either pi pulse on the nuclei or a delay of the same length for reference
    rho_a=shaped_pulse_af(spin_system,L,Nx,Ny,rho0,parameters.n_frq(n),parameters.n_pwr,...
                                                   parameters.n_dur,   parameters.n_phi,...
                                                   parameters.n_rnk,   parameters.method); %#ok<PFBNS>
    rho_b=shaped_pulse_af(spin_system,L,Nx,Ny,rho0,parameters.n_frq(n),0*parameters.n_pwr,...
                                                   parameters.n_dur,     parameters.n_phi,...
                                                   parameters.n_rnk,     parameters.method);

    % Pi/2 pulse on the electron
    rho_a=shaped_pulse_af(spin_system,L,Ex,Ey,rho_a,parameters.e_frq,  parameters.e_pwr,...
                                                    parameters.e_dur/2,parameters.e_phi,...
                                                    parameters.e_rnk,  parameters.method);
    rho_b=shaped_pulse_af(spin_system,L,Ex,Ey,rho_b,parameters.e_frq,  parameters.e_pwr,...
                                                    parameters.e_dur/2,parameters.e_phi,...
                                                    parameters.e_rnk,  parameters.method);

    % Optionally, run the spin echo stage
    if isfield(parameters,'tau')&&(parameters.tau~=0)
        
        % First evolution period
        rho_a=evolution(spin_system,L,[],rho_a,parameters.tau,1,'final');
        rho_b=evolution(spin_system,L,[],rho_b,parameters.tau,1,'final');
    
        % Pi pulse on the electron
        rho_a=shaped_pulse_af(spin_system,L,Ex,Ey,rho_a,parameters.e_frq,parameters.e_pwr,...
                                                        parameters.e_dur,parameters.e_phi,...
                                                        parameters.e_rnk,parameters.method);
        rho_b=shaped_pulse_af(spin_system,L,Ex,Ey,rho_b,parameters.e_frq,parameters.e_pwr,...
                                                        parameters.e_dur,parameters.e_phi,...
                                                        parameters.e_rnk,parameters.method);  
                                                            
        % Second evolution period
        rho_a=evolution(spin_system,L,[],rho_a,parameters.tau,1,'final');
        rho_b=evolution(spin_system,L,[],rho_b,parameters.tau,1,'final');
        
    end
    
    % Return the ratio of RF-on and RF-off
    answer(n)=(parameters.coil'*rho_a)/(parameters.coil'*rho_b);
                                           
end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available in Liouville space.');
end
if ~isfield(parameters,'method')
    error('shaped pulse simulation method must be specified in parameters.method field.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'e_frq')
    error('electron pulse frequency must be specified in parameters.e_frq field.');
end
if (~isnumeric(parameters.e_frq))||(~isreal(parameters.e_frq))||(~isscalar(parameters.e_frq))
    error('parameters.e_frq must be a real scalar.');
end
if ~isfield(parameters,'e_phi')
    error('electron pulse phase must be specified in parameters.e_phi field.');
end
if (~isnumeric(parameters.e_phi))||(~isreal(parameters.e_phi))||(~isscalar(parameters.e_phi))
    error('parameters.e_phi must be a real scalar.');
end
if ~isfield(parameters,'e_pwr')
    error('electron pulse power must be specified in parameters.e_pwr field.');
end
if (~isnumeric(parameters.e_pwr))||(~isreal(parameters.e_pwr))||...
   (~isscalar(parameters.e_pwr))||(parameters.e_pwr<=0)
    error('parameters.e_pwr must be a positive real scalar.');
end
if ~isfield(parameters,'e_dur')
    error('electron pulse duration must be specified in parameters.e_dur field.');
end
if (~isnumeric(parameters.e_dur))||(~isreal(parameters.e_dur))||...
   (~isscalar(parameters.e_dur))||(parameters.e_dur<=0)
    error('parameters.e_dur must be a positive real scalar.');
end
if ~isfield(parameters,'e_rnk')
    error('electron pulse grid rank must be specified in parameters.e_rnk field.');
end
if (~isnumeric(parameters.e_rnk))||(~isreal(parameters.e_rnk))||...
   (~isscalar(parameters.e_rnk))||(mod(parameters.e_rnk,1)~=0)||...
   (parameters.e_rnk<1)
    error('parameters.e_rnk must be a positive real integer.');
end
if ~isfield(parameters,'n_frq')
    error('nuclear pulse frequencies must be specified in parameters.n_frq field.');
end
if (~isnumeric(parameters.n_frq))||(~isreal(parameters.n_frq))
    error('parameters.n_frq must be an array of real numbers.');
end
if ~isfield(parameters,'n_phi')
    error('nuclear pulse phase must be specified in parameters.n_phi field.');
end
if (~isnumeric(parameters.n_phi))||(~isreal(parameters.n_phi))||(~isscalar(parameters.n_phi))
    error('parameters.n_phi must be a real scalar.');
end
if ~isfield(parameters,'n_pwr')
    error('nuclear pulse power must be specified in parameters.n_pwr field.');
end
if (~isnumeric(parameters.n_pwr))||(~isreal(parameters.n_pwr))||...
   (~isscalar(parameters.n_pwr))||(parameters.n_pwr<=0)
    error('parameters.n_pwr must be a positive real scalar.');
end
if ~isfield(parameters,'n_dur')
    error('nuclear pulse duration must be specified in parameters.n_dur field.');
end
if (~isnumeric(parameters.n_dur))||(~isreal(parameters.n_dur))||...
   (~isscalar(parameters.n_dur))||(parameters.n_dur<=0)
    error('parameters.n_dur must be a positive real scalar.');
end
if ~isfield(parameters,'n_rnk')
    error('nuclear pulse grid rank must be specified in parameters.n_rnk field.');
end
if (~isnumeric(parameters.n_rnk))||(~isreal(parameters.n_rnk))||...
   (~isscalar(parameters.n_rnk))||(mod(parameters.n_rnk,1)~=0)||...
   (parameters.n_rnk<1)
    error('parameters.n_rnk must be a positive real integer.');
end
end

% "English is hard, but we can it!"
%
% A poster seen at the Chemistry Department
% of Munich University of Technology

