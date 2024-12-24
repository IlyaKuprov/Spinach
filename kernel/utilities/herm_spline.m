% Cubic Hermite spline on [0,1] interval from values and deriva-
% tives at the interval edges. Syntax:
%
%                  y=herm_spline(f0,df0,f1,df1,x)
%
% Parameters:
%
%    f0  - function value(s) at the left edge, a real
%          scalar or array
%
%    df0 - function derivative(s) at the left edge,
%          a real scalar or array
%
%    f1  - function value(s) at the right edge, a real
%          scalar or array
%
%    df1 - function derivative(s) at the right edge,
%          a real scalar or array
%
%    x   - query point(s) inside [0,1] interval,
%          a real scalar or array
%
% Function values and derivatives can be scalars (in which case
% the same spline is evaluated at all query points) or arrays of
% the same size as x, in which case multiple splines are evalua-
% ted at their corresponding query points.
%
% Outputs:
%
%    y   - the value of the spline(s) at the query point(s)
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=herm_spline.m>

function y=herm_spline(f0,df0,f1,df1,x)

% Check consistency
grumble(f0,df0,f1,df1,x);

% Adapt to input
if isscalar(f0)&&isscalar(df0)&&...
   isscalar(f1)&&isscalar(df1)&&(~isscalar(x))
    f0=f0*ones(size(x)); df0=df0*ones(size(x));
    f1=f1*ones(size(x)); df1=df1*ones(size(x));
end

% Preallocate the output
y=zeros(size(x));

% Loop over the entries
for n=1:numel(x)

    % Get spline coefficients (x^3 -> x^0)
    c=[ 1   2   1  -2;
       -2  -3  -1   3;
        1   0   0   0;
        0   1   0   0]*[df0(n); f0(n); df1(n); f1(n)];
           
    % Evaluate at the query point
    y(n)=c(1)*x(n)^3+c(2)*x(n)^2+c(3)*x(n)+c(4);

end

end

% Consistency enforcement
function grumble(f0,df0,f1,df1,x)
if (~isnumeric(f0))||(~isreal(f0))||any(~isfinite(f0),'all')||...
   (~isnumeric(df0))||(~isreal(df0))||any(~isfinite(df0),'all')||...
   (~isnumeric(f1))||(~isreal(f1))||any(~isfinite(f1),'all')||...
   (~isnumeric(df1))||(~isreal(df1))||any(~isfinite(df1),'all')||...
   (~isnumeric(x))||(~isreal(x))||any(~isfinite(x),'all')
    error('all inputs must be numeric, real, and finite.');
end
if ((~isscalar(f0))&&(~all(size(f0)==size(x))))||...
   ((~isscalar(df0))&&(~all(size(df0)==size(x))))||...
   ((~isscalar(f1))&&(~all(size(f1)==size(x))))||...
   ((~isscalar(df1))&&(~all(size(df1)==size(x))))
    error('f0,f1,df0,df1 arrays must either be scalars, or have the same size as x.');
end
end

% I was playing in a tournament in Germany one year when a man approached
% me. Thinking he just wanted an autograph, I reached for my pen, when the
% man made a startling announcement... "I've solved chess!" I sensibly
% started to back away in case the man was dangerous as well as insane, but
% the man continued: "I'll bet you 50 marks that if you come back to my
% hotel room I can prove it to you." Well, 50 marks was 50 marks, so I
% humored the fellow and accompanied him to his room. Back at the room, we
% sat down at his chess board. "I've worked it all out, white mates in 12
% moves no matter what." I played with black perhaps a bit incautiously,
% but I found to my horror that white's pieces coordinated very strangely,
% and that I was going to be mated on the 12th move! I tried again, and I
% played a completely different opening that couldn't possibly result in
% such a position, but after a series of very queer-looking moves, once
% again I found my king surrounded, with mate to fall on the 12th move. I
% asked the man to wait while I ran downstairs and fetched Emmanuel Lasker,
% who was world champion before me. He was extremely skeptical, but agreed
% to at least come and play. Along the way we snagged Alekhine, who was
% then world champion, and the three of us ran back up to the room.
% 
% Lasker took no chances, but played as cautiously as could be, yet after
% a bizarre, pointless-looking series of maneuvers, found himself hemmed
% in a mating net from which there was no escape. Alekhine tried his hand,
% too, but all to no avail.
% 
% It was awful! Here we were, the finest players in the world, men who had
% devoted our very lives to the game, and it was all over! The tournaments,
% the matches, everything - chess had been solved, white wins.
% 
% We killed him, of course.
% 
% Jose Raul Capablanca, only half joking.

