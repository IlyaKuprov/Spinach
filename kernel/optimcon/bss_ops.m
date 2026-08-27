% Bloch-Siegert response operators for the optimal control module. For
% each control channel, returns the operator whose coefficient in every
% time slice of a GRAPE optimisation is the square of the physical con-
% trol amplitude on that channel. The operator collects the second-order
% Bloch-Siegert frequency shifts of every spin in the system:
%
%          B=sum_n (gamma_n/gamma_c)^2*[1/(2*(omega_n+omega_c))
%                +(foreign isotopes only) 1/(2*(omega_n-omega_c))]*Lz_n
%
% where omega_n are signed laboratory frame Zeeman frequencies from
% spin_system.inter.basefrqs and omega_c is the signed carrier frequen-
% cy of the channel. Spins belonging to the isotope that the channel
% addresses receive only the never-resonant term because their resonant
% term is the control operator itself, which GRAPE propagates exactly.
% Foreign isotopes receive both terms, which sum into the Ramsey shift
% omega_n*omega_1n^2/(omega_n^2-omega_c^2). Syntax:
%
%       resp_ops=bss_ops(spin_system,channels,carrier_frq)
%
% Parameters:
%
%    channels     - cell array of isotope strings, one per
%                   control operator, e.g. {'1H','1H'} for
%                   an X,Y control operator pair
%
%    carrier_frq  - row vector of signed carrier frequencies
%                   (rad/s), one per control operator; for a
%                   transmitter on resonance this is the value
%                   of inter.basefrqs for the channel isotope
%
% Outputs:
%
%    resp_ops     - cell array of Bloch-Siegert response ope-
%                   rators, one per control operator, in the
%                   formalism of the spin system provided
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bss_ops.m>

function resp_ops=bss_ops(spin_system,channels,carrier_frq)

% Check consistency
grumble(spin_system,channels,carrier_frq);

% Preallocate the answer
resp_ops=cell(size(channels));

% Loop over control channels
for c=1:numel(channels)

    % Get the channel magnetogyric ratio
    gamma_chan=spin(channels{c});

    % Preallocate the response operator
    resp_op=mprealloc(spin_system,1);

    % Loop over the spins
    for n=1:spin_system.comp.nspins

        % Skip non-spin particles
        if ~strcmp(spin_system.comp.types{n},'S'), continue; end

        % Field scaling factor and Zeeman frequency
        scal=(spin_system.inter.gammas(n)/gamma_chan)^2;
        omega_n=spin_system.inter.basefrqs(n);

        % Never-resonant term for every spin
        coeff=scal/(2*(omega_n+carrier_frq(c)));

        % Resonant-capable term for foreign isotopes
        if ~strcmp(spin_system.comp.isotopes{n},channels{c})
            coeff=coeff+scal/(2*(omega_n-carrier_frq(c)));
        end

        % Add the weighted longitudinal operator
        resp_op=resp_op+coeff*operator(spin_system,{'Lz'},{n});

    end

    % Store without clean-up, elements scale as inverse Zeeman frequencies
    resp_ops{c}=resp_op;

end

end

% Consistency enforcement
function grumble(spin_system,channels,carrier_frq)
if ~isfield(spin_system,'comp')
    error('spin_system must be a Spinach data structure.');
end
if (~iscell(channels))||any(~cellfun(@ischar,channels(:)))
    error('channels must be a cell array of character strings.');
end
if (~isnumeric(carrier_frq))||(~isreal(carrier_frq))||...
   (~all(isfinite(carrier_frq)))||(numel(carrier_frq)~=numel(channels))
    error('carrier_frq must be a finite real vector with one element per channel.');
end
for c=1:numel(channels)
    if ~ismember(channels{c},spin_system.comp.isotopes)
        error(['channel isotope ' channels{c} ' is not present in the system.']);
    end
    if carrier_frq(c)==0
        error('zero carrier frequency, Bloch-Siegert denominators are undefined.');
    end
    for n=1:spin_system.comp.nspins
        if ~strcmp(spin_system.comp.types{n},'S'), continue; end
        omega_n=spin_system.inter.basefrqs(n);
        if abs(omega_n+carrier_frq(c))<1e-6*abs(carrier_frq(c))
            error('degenerate omega_n+omega_c denominator in the system.');
        end
        if (~strcmp(spin_system.comp.isotopes{n},channels{c}))&&...
           (abs(omega_n-carrier_frq(c))<1e-6*abs(carrier_frq(c)))
            error('degenerate omega_n-omega_c denominator in the system.');
        end
    end
end
end

% There is nothing noble in being superior to your fellow
% man; true nobility is being superior to your former self.
%
% Ernest Hemingway

