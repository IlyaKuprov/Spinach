% Three-pulse DEER pulse sequence. Idealized hard pulses are used,
% each pulse only affects its specific electron or transition, de-
% pending on the pulse operators supplied. Syntax:
%
%        deer=deer_3p_hard_deer(spin_system,parameters,H,R,K)
%
% Parameters:
%
%   parameters.rho0            initial state
%
%   parameters.coil_prob       detection state on probe spin
%
%   parameters.stepsize        increment time for the pump pulse
%                              sandwich
%
%   parameters.nsteps          number of steps for the pump pulse
%                              sandwich
%
%   parameters.ex_prob         excitation operators to be used for
%   parameters.ex_pump         the probe and pump electron respec-
%                              tively.
%
%   parameters.output          'brief' returns just the DEER trace,
%                              'detailed' also returns excitation 
%                               profiles and the EPR spectrum.
%
%   H  - Hamiltonian matrix, received from context function
%
%   R  - relaxation superoperator, received from context function
%
%   K  - kinetics superoperator, received from context function
%
% If 'detailed' is selected as the output option, the following pa-
% rameters are also required:
%
%   parameters.ex_hard          hard pulse excitation operator
%
%   parameters.spectrum_sweep   sweep width of the EPR spectrum, Hz
%
%   parameters.spectrum_nsteps  number of time steps in the FID
%
%   parameters.coil_pump        detection state on pump spin
%
% Outputs:
%
%   deer.hard_pulse_fid  - ('detailed') free induction decay
%                          after a non-selective ideal pulse
%
%   deer.prob_pulse_fid  - ('detailed') free induction decay
%                          after just the the probe pulse
%
%   deer.pump_pulse_fid  - ('detailed') free induction decay
%                          after just the pump pulse
%
%   deer.deer_trace      - DEER signal
%
% Note: hard pulses are only appropriate for spin-1/2 systems; for
%       higher spin systems transition selective pulse operators
%       must be supplied.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
% nurit.manukovsky@weizmann.ac.il
% daniella.goldfarb@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=deer_3p_hard_deer.m>

function deer=deer_3p_hard_deer(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Compute detailed output if required
if strcmp(parameters.output,'detailed')
    
    % Apply a hard 90-degree pulse
    rho_hard=step(spin_system,parameters.ex_hard,parameters.rho0,pi/2);
    
    % Return the free induction decay
    deer.hard_pulse_fid=evolution(spin_system,L,parameters.coil_pump+parameters.coil_prob,rho_hard,...
                        1/parameters.spectrum_sweep,parameters.spectrum_nsteps,'observable');
    
    % Apply a selective 90-degree pulse on the probe spin
    rho=step(spin_system,parameters.ex_prob,parameters.rho0,pi/2);
    
    % Return the free induction decay on the probe spin
    deer.prob_pulse_fid=evolution(spin_system,L,parameters.coil_prob,rho,1/parameters.spectrum_sweep,...
                                  parameters.spectrum_nsteps,'observable');
    
    % Apply a selective 180-degree pulse on the pump spin
    rho=step(spin_system,parameters.ex_pump,parameters.rho0,pi/2);
    
    % Return the free induction decay on the pump spin
    deer.pump_pulse_fid=evolution(spin_system,L,parameters.coil_pump,rho,1/parameters.spectrum_sweep,...
                                  parameters.spectrum_nsteps,'observable');
    
end

% First pulse
rho=step(spin_system,parameters.ex_prob,parameters.rho0,pi/2);

% Evolution
rho_stack=evolution(spin_system,L,[],rho,parameters.stepsize,parameters.nsteps,'trajectory');

% Second pulse
rho_stack=step(spin_system,parameters.ex_pump,rho_stack,pi);
     
% Evolution
rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),parameters.stepsize,parameters.nsteps,'refocus');

% Third pulse
rho_stack=step(spin_system,parameters.ex_prob,rho_stack,pi);

% Evolution
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.stepsize,parameters.nsteps,'final');

% Observation
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    deer.deer_trace=cellfun(@(x)(trace(parameters.coil_prob'*x)),rho_stack)/norm(parameters.coil_prob,'fro');
else
    deer.deer_trace=parameters.coil_prob'*rho_stack/norm(parameters.coil_prob,2);
end

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'rho0')
    error('the initial state must be provided in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil_prob')
    error('the detection state on the probe spin must be provided in parameters.coil_prob variable.');
end
if ~isfield(parameters,'stepsize')
    error('increment time for the pump pulse sandwich must be provided in parameters.stepsize variable.');
end
if (~isnumeric(parameters.stepsize))||(~isreal(parameters.stepsize))||...
   (numel(parameters.stepsize)~=1)||(parameters.stepsize<=0)
    error('parameters.stepsize must be a positive real number.');
end
if ~isfield(parameters,'nsteps')
    error('number of steps for the pump pulse sandwich must be provided in parameters.nsteps variable.');
end
if (~isnumeric(parameters.nsteps))||(~isreal(parameters.nsteps))||...
   (numel(parameters.nsteps)~=1)||(mod(parameters.nsteps,1)~=0)||(parameters.nsteps<1)
    error('parameters.nsteps must be a positive real integer.');
end
if ~isfield(parameters,'ex_prob')
    error('probe pulse operator must be provided in parameters.ex_prob variable.');
end
if ~isfield(parameters,'ex_pump')
    error('pump pulse operator must be provided in parameters.ex_pump variable.');
end
if ~isfield(parameters,'output')
    error('output type must be provided in parameters.output variable.');
end
if ~ischar(parameters.output)
    error('parameters.output must be a character string.');
end
if ~ismember(parameters.output,{'brief','detailed'})
    error('the available values for parameters.output are ''brief'' and ''detailed''.');
end
if strcmp(parameters.output,'detailed')
    if ~isfield(parameters,'ex_hard')
        error('hard pulse excitation operator must be provided in parameters.ex_hard for detailed output.');
    end
    if ~isfield(parameters,'spectrum_sweep')
        error('sweep width of the EPR spectrum must be provided in parameters.spectrum_sweep for detailed output.');
    end
    if (~isnumeric(parameters.spectrum_sweep))||(~isreal(parameters.spectrum_sweep))||...
       (numel(parameters.spectrum_sweep)~=1)||(parameters.spectrum_sweep<=0)
        error('parameters.spectrum_sweep must be a positive real number.');
    end
    if ~isfield(parameters,'spectrum_nsteps')
        error('number of time steps in the fid must be provided in parameters.spectrum_nsteps for detailed output.');
    end
    if (~isnumeric(parameters.spectrum_nsteps))||(~isreal(parameters.spectrum_nsteps))||...
       (numel(parameters.spectrum_nsteps)~=1)||(mod(parameters.spectrum_nsteps,1)~=0)||...
       (parameters.spectrum_nsteps<1)
        error('parameters.spectrum_nsteps must be a positive real integer.');
    end
    if ~isfield(parameters,'coil_pump')
        error('detection state on pump spin must be provided in parameters.coil_pump for detailed output.');
    end
end
end

% On two occasions I have been asked [by members of Parliament], "Pray,
% Mr. Babbage, if you put into the machine wrong figures, will the right
% answers come out?" I am not able rightly to apprehend the kind of
% confusion of ideas that could provoke such a question.
%
% Charles Babbage

