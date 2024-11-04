% MD5 hash of any Matlab object as a hex string. Identical sparse 
% and full matrices return different hashes. Syntax:
%
%                       hashstr=md5_hash(A)
%
% Parameters:
%
%     A            - Matlab object of any type
%
% Outputs:
%
%     hashstr      - hexadecimal character string  
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=md5_hash.m>

function hashstr=md5_hash(A)

% Make a bytestream
A=serializeToBytes(A);

% Compute MD5 hash
hashstr=digestMD5(A);

% Convert into a hex string
hashstr=sprintf('%.2x',hashstr);

end

% The basic principle of the new education is to be that dunces and
% idlers must not be made to feel inferior to intelligent and indus-
% trious pupils. That would be "undemocratic". These differences be-
% tween pupils - for there are obviously and nakedly individual dif-
% ferences - must be disguised. This can be done at various levels.
% At universities, examinations must be framed so that nearly all the
% students get good marks. [...] But all the time there must be no
% faintest hint that they are inferior to the children who are at 
% work. Whatever nonsense they are engaged in must have - I believe
% the English already use the phrase - "parity of esteem".
%
% C.S. Lewis

