% Converts the Josephson and charging energies of a transmon into
% the Duffing oscillator frequency and anharmonicity expected by
% the bosonic mode specification interface of create.m using the
% asymptotic transmon expressions (Koch et al., https://doi.org/
% 10.1103/PhysRevA.76.042319):
%
%          frq=sqrt(8*ej*ec)-ec,     anharm=-ec
%
% Syntax:
%
%               [frq,anharm]=ejec2duffing(ej,ec)
%
% Parameters:
%
%   ej      - Josephson energies in Hz (energy over the
%             Planck constant), an array of positive
%             real numbers
%
%   ec      - charging energies in Hz (energy over the
%             Planck constant), an array of positive
%             real numbers of the same size as ej
%
% Outputs:
%
%   frq     - transition frequencies in Hz, to be placed
%             into inter.modes.frqs
%
%   anharm  - Duffing anharmonicities in Hz, to be placed
%             into inter.modes.anharms
%
% Note: the asymptotic expressions are only accurate deep in the
%       transmon regime ej/ec>>1; a warning is issued when the
%       ratio is smaller than 20.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ejec2duffing.m>

function [frq,anharm]=ejec2duffing(ej,ec)

% Check consistency
grumble(ej,ec);

% Warn outside the transmon regime
if any(ej./ec<20,'all')
    warning('ej/ec ratio below 20, asymptotic transmon expressions are inaccurate.');
end

% Run the conversion
frq=sqrt(8*ej.*ec)-ec; anharm=-ec;

end

% Consistency enforcement
function grumble(ej,ec)
if (~isnumeric(ej))||(~isreal(ej))||any(~isfinite(ej),'all')||any(ej<=0,'all')
    error('ej must be an array of positive real numbers.');
end
if (~isnumeric(ec))||(~isreal(ec))||any(~isfinite(ec),'all')||any(ec<=0,'all')
    error('ec must be an array of positive real numbers.');
end
if ~isequal(size(ej),size(ec))
    error('ej and ec must have the same size.');
end
end

% For a successful technology, reality must take precedence
% over public relations, for Nature cannot be fooled.
%
% Richard Feynman

