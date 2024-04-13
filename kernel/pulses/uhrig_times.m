% Uhrig's UDD decoupling sequence timings. Syntax:
%
%               time_delays=uhrig_times(T,N)
%
% Parameters:
%
%    T - total duration of the sequence (sum of all delays)
%
%    N - number of pulses in the sequence
%
% Outputs:
%
%    time_delays - list of delays between ideal pulses in
%                  the UDD sequence; the first pulse goes 
%                  after the first delay, and there is a 
%                  delay after the last pulse
%
% i.kuprov@soton.ac.uk
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=uhrig_times.m>

function time_delays=uhrig_times(T,N)

% Check consistency
grumble(T,N);

% Use the formula from WSW's 2009 JCP paper
time_positions=T*(sin(pi*(1:N)/(2*N+2)).^2-0.5);

% Convert positions to delays
time_delays=diff(time_positions);

% Add the starting and the trailing delay
chunk=(T-sum(time_delays))/2;
time_delays=[chunk time_delays chunk];

end

% Consistency enforcement
function grumble(T,N)
if (~isreal(T))||(~isscalar(T))||(T<=0)
    error('T must be a positive real scalar.');
end
if (~isreal(N))||(~isscalar(N))||(N<1)||(mod(N,1)~=0)
    error('N must be a positve real integer.');
end
end

% "In a legislative effort to force vendors to hand over the plaintext
% contents of encrypted communications, Prime Minister Malcom Turnbull 
% was confronted with the pesky laws of mathematics; that is, that with
% properly implemented crypto, no one, including vendors, could simply
% decrypt the user data. Turnbull replied 'Well the laws of Australia
% prevail in Australia, I can assure you of that. The laws of mathema-
% tics are very commendable, but the only law that applies in Australia
% is the law of Australia'. Finally, some progress in the war on math!"
%
% Laudation for Malcolm Turnbull's Pwnie Award nomination

