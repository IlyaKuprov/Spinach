% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. Syntax:
%
%                       kcolourbar(x)
%
% Parameters:
%
%    x - a character string
%
% Outputs:
% 
%    creates or updates the colour bar 
%    in the current axis system
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kcolourbar.m>

function kcolourbar(x)

% Check consistency
grumble(x);

% Ticks to LaTeX
cb=colorbar('TickLabelInterpreter','latex',...
            'FontSize',12);

% Label to LaTex
cb.Label.Interpreter='latex';
cb.Label.FontSize=13;
cb.Label.String=x; 

end

% Consistency enforcement
function grumble(x)
if ~ischar(x)
    error('x must be a character string');
end
end

% I've often remarked that identity politics is the product of 
% prosperity. The movement's nitpicking about undetectable-to-
% the-human-eye 'microaggressions' is an indulgence, most com-
% monly among affluent white people hungry for the illusion of
% having problems (and enemies) for a sense of meaning. [...]
%
% Arguably, humanity by nature requires challenging obstacles 
% to overcome in order to feel purposeful, so that when socie-
% ties become too coddling, too safe and too comfortable, they
% also become neurotic: they construct make-believe problems
% to chafe against, like scratching-posts for pets.
%
% Lionel Shriver

