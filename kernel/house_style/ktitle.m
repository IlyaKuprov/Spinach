% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. Syntax:
%
%                          ktitle(x)
%
% Parameters:
%
%    x - a character string
%
% Outputs:
% 
%    creates or updates the title in 
%    the current axis system
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ktitle.m>

function ktitle(x)

% Check consistency
grumble(x);

% Bold title rendered by LaTeX
title(['\textbf{' x '}'],'Interpreter','latex');

end

% Consistency enforcement
function grumble(x)
if ~ischar(x)
    error('x must be a character string');
end
end

% 3:58 Признав иных, я вслед за тем в одном
%      Узнал того, кто от великой доли
%      Отрекся в малодушии своем.
% 
% 3:61 И понял я, что здесь вопят от боли
%      Ничтожные, которых не возьмут
%      Ни бог, ни супостаты божьей воли.
%
% Данте Алигьери, Божественная Kомедия,
% перевод М. Лозинского

