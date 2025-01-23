% R sequence compiler. Uses the fact that R-sequences are very
% repetitive to pre-compile the minimal number of pulse propa-
% gators. Syntax:
%
%   [P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,...
%                       pulse_amp,pulse_dur,element_type)
%
% Parameters:
%
%    L         - background Liouvillian
%
%    Sx,Sy     - Cartesian spin operators pertaining to
%                the spins affected by the pulses
%
%    pulse_phi - the sequence of pulse phases, radians
%
%    pulse_amp - RF nutation frequency in rad/s, a scalar 
%                because R-sequences are phase-modulated
%
%    pulse_dur - duration of the pulses in the sequence 
%                element, a vector with the length mat-
%                ching the number of pulses in the sequ-
%                ence element (seconds)
%
%    element_type - R element needs to be an inversion 
%                   pulse; common ones are:
%
%                    '180_pulse'   : simple inversion pulse
%
%                    '90270_pulse' : composite inversion pulse
%
% Outputs:
%
%    P - unique propagators, a cell array of matrices
%
%    T - an index array of the same dimension as pulse_phi,
%        specifying which propagator is to be used at which
%        slice of the phase sequence
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rseq_compiler.m>
 
function [P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,...
                             pulse_amp,pulse_dur,element_type)

% Check consistency
grumble(L,Sx,Sy,pulse_phi,pulse_amp,pulse_dur,element_type);

% Find unique phases
[phi,~,T]=unique(pulse_phi);

% Preallocate propagators
P=cell(numel(phi),1);

% Compute propagators
switch element_type
    
    % Sequence of pi pulses
    case '180_pulse'
        
        % Run matrix exponentials
        for n=1:numel(phi)
            LP=L+pulse_amp*(Sx*cos(phi(n))+Sy*sin(phi(n)));
            P{n}=propagator(spin_system,LP,pulse_dur);
        end
    
    % Sequence of composite pulses
    case '90270_pulse'
        
        % Run matrix exponentials
        for n=1:2:numel(phi)
            LP=L+pulse_amp*(Sx*cos(phi(n))+Sy*sin(phi(n)));
            P{n}=propagator(spin_system,LP,pulse_dur(1));
            LP=L+pulse_amp*(Sx*cos(phi(n+1))+Sy*sin(phi(n+1)));
            P{n+1}=propagator(spin_system,LP,pulse_dur(2));
        end
        
    otherwise
        
        % Complain and bomb out
        error('Unknown pulse element type.');
        
end

end

% Consistency enforcement
function grumble(L,Sx,Sy,pulse_phi,pulse_amp,pulse_dur,element_type)
if ~ischar(element_type)
    error('element_type must be a character string.');
end
if (~isnumeric(L))||(~isnumeric(Sx))||(~isnumeric(Sy))
    error('L, Sx, Sy must be matrices.');
end
if (~ishermitian(Sx))||(~ishermitian(Sy))
    error('Sx and Sy matrices must be Hermitian.');
end
if (~isnumeric(pulse_amp))||(~isreal(pulse_amp))||...
   (~isscalar(pulse_amp))||(pulse_amp<0)
    error('pulse_amp must be a non-negative real scalar.');
end
if (~isnumeric(pulse_phi))||(~isreal(pulse_phi))
    error('pulse_phi must be a real numeric array.');
end
if (~isnumeric(pulse_dur))||(~isreal(pulse_dur))
    error('pulse_dur must be a real numeric array.');
end
end

% The more I learn about people, the more I like my dog.
%
% Mark Twain

