% R-sequences described in Malcolm Levitt's review:
%
%          https://doi.org/10.1002/9780470034590.emrstm0551
%
% Nomenclature is based on the following notation RN_{n}^{\nu}. Syntax:
%
%  [phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,...
%                                         phase_factor,n_cycle_repeats,mas_rate,...
%                                         element_type,supercycle_type)
%
% Parameters:
%
%   n_rotor_periods     -  "small n" symmetry number, gives number of
%                           rotor periods required in the R symmetry
%
%   n_blocks_per_period -  "capital n" symmetry number, gives the number of
%                           R elements contained within the R symmetry
%
%   phase_factor        -  "nu" to calculate the alternating phase in the 
%                           R sequence:
%
%                          180*nu/N = 180*phase_factor/n_blocks_per_period
%
%   n_cycle_repeats     -  number of times the full R sequence is applied 
%                 
%   mas_rate            -  rotor spinning rate, Hz
%
%   element_type        -  R element needs to be an inversion pulse. Common
%                          R elements are:
%
%                           '180_pulse'   : simple inversion pulse
%
%                           '90270_pulse' : composite inversion pulse
%
%   supercycle_type     -  The R sequence can be repeated multiple
%                          times in combination with supercycles, for
%                          improved performance, removal of undesired 
%                          higher order terms. If the he unmodified R
%                          is denoted [phase], this can either be in-
%                          verted, [-phase], or have an overall phase
%                          added to it, [phase]_addph.
%
%                          Common supercycles are:
%
%                           'hetero_single_quantum' 
%                           [phase]_0:[-phase]_0:[phase]_120:[-phase]_120:[phase]_240:[-phase]_240
%
%                           'homo_double_quantum_nucycle'
%                           [phase]_0:[-phase]_0
%
%                           'homo_double_quantum_nupicycle'
%                           [phase]_0:[-phase]_0:[-phase]_180:[phase]_180
%
% Output:
%  
%   phases    - the sequence of pulse phases, radians
%
%   pulse_amp - RF nutation frequency in rad/s, a scalar because
%               R-sequences are phase-modulated
%
%   pulse_dur - duration of the pulses in the sequence element,
%               a vector with the length matching the number of
%               pulses in the sequence element (seconds)
%               
% i.kuprov@soton.ac.uk
% m.carravetta@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rsequence.m>

function [phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,...
                                                phase_factor,n_cycle_repeats,mas_rate,...
                                                element_type,supercycle_type)
% Check consistency
grumble(n_rotor_periods,n_blocks_per_period,phase_factor,...
        n_cycle_repeats,mas_rate,element_type,supercycle_type);

% Get rotor period (seconds)
rotor_period=1/mas_rate;

% Get R element (either simple or composite) duration (seconds)
r_element_dur=n_rotor_periods*rotor_period/n_blocks_per_period;

% Decide sequence element
switch element_type

    % A train of pi pulses
    case '180_pulse'

        % Get RF amplitude (rad/s)
        pulse_amp=pi/r_element_dur;

        % One duration (seconds)
        pulse_dur=r_element_dur;

        % Calculate phase grid
        phases=zeros(n_blocks_per_period,1);
        for q=0:(n_blocks_per_period-1)
            phases(q+1)=((-1)^q)*pi*phase_factor/n_blocks_per_period;
        end
    
    % A train of composite pulses with 
    % the overall duration of a 360 pulse
    case '90270_pulse'

        % Get RF amplitude (rad/s)
        pulse_amp=2*pi/r_element_dur;

        % Two durations (seconds)
        pulse_dur(1)=(1/4)*r_element_dur;
        pulse_dur(2)=(3/4)*r_element_dur;

        % Calculate phase grid
        phases=zeros(n_blocks_per_period*2,1);
        for q=0:(n_blocks_per_period-1)
            phases(2*q+1)=((-1)^q)*pi*phase_factor/n_blocks_per_period;
            phases(2*q+2)=phases(2*q+1)+pi;
        end

    otherwise

        % Complain and bomb out
        error('unknown sequence element type.');

end

% Apply the supercycle
switch supercycle_type

    case 'hetero_single_quantum'

        % Replicate with inverted phases
        phases=[phases; -phases];

        % Replicate with shifted phases
        phases=[phases; phases+2*pi/3; phases+4*pi/3];

    case 'homo_double_quantum_nupicycle'

        % Replicate with inverted phases 
        % and a two-step pi supercycle
        phases=[phases; -phases; -phases+pi; phases+pi];

    case 'homo_double_quantum_nucycle'

        % Replicate with inverted phases 
        phases=[phases; -phases];

    otherwise

        % Complain and bomb out
        error('unknown supercycle type.');

end
   
% Repeat the specified number of times
phases=kron(ones(n_cycle_repeats,1),phases);

end

% Consistency enforcement
function grumble(n_rotor_periods,n_blocks_per_period,phase_factor,...
                 n_cycle_repeats,mas_rate,element_type,supercycle_type)
if (~isnumeric(n_rotor_periods))||(~isreal(n_rotor_periods))||...
   (~isscalar(n_rotor_periods))||(mod(n_rotor_periods,1)~=0)||...
   (n_rotor_periods<1)
    error('n_rotor_periods must be a positive real scalar.');
end
if (~isnumeric(n_blocks_per_period))||(~isreal(n_blocks_per_period))||...
   (~isscalar(n_blocks_per_period))||(mod(n_blocks_per_period,1)~=0)||...
   (n_blocks_per_period<1)
    error('n_blocks_per_period must be a positive real scalar.');
end
if (~isnumeric(n_cycle_repeats))||(~isreal(n_cycle_repeats))||...
   (~isscalar(n_cycle_repeats))||(mod(n_cycle_repeats,1)~=0)||...
   (n_cycle_repeats<1)
    error('n_cycle_repeats must be a positive real scalar.');
end
if (~isnumeric(mas_rate))||(~isreal(mas_rate))||(~isscalar(mas_rate))
    error('mas_rate must be a real scalar.');
end
if (~isnumeric(phase_factor))||(~isreal(phase_factor))||...
   (~isscalar(phase_factor))
    error('phase_factor must be a real scalar.');
end
if ~ischar(element_type)
    error('element_type must be a character string.');
end
if ~ischar(supercycle_type)
    error('supercycle_type must be a character string.');
end
end

% Effortlessness is hard.
%
% Sarah Sands

