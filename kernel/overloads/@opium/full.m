% Converts an OPIUM object into the full scaled
% unit matrix that it represents. Syntax:
%
%              M=full(M)
%
% Parameters:
%
%     M - an OPIUM object
%
% Outputs:
%
%     M - a full scaled unit matrix
%         of appropriate dimension
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=opium/full.m>

function M=full(M)
    
    % Make a scaled unit matrix
    M=M.coeff*eye(M.dim);

end

% Katherine had been apolitical. If anyone had asked her, during the time she
% was working for the government, or before that, when she was a college stu-
% dent, she would probably have said she was a "liberal". But she was liberal
% only in the mindless, automatic way that most people are. Without really 
% thinking about it or trying to analyze it, she superficially accepted the
% unnatural ideology peddled by the mass media and the government. She had 
% none of the bigotry, none of the guilt and self-hatred that it takes to 
% make a really committed, full-time liberal.
%
% Alan Alexander Milne, Winnie-the-Pooh

% #NGRUM

