% Performs arbitrary-order finite-difference differentiation of a
% user-supplied row or column vector. Uses central finite-differe-
% nce stencils in the middle and sided stencils of the same order 
% of accuracy on the sides. Syntax:
%
%                     dx=fdvec(x,npoints,order)
%
% Parameters:
%
%    x       - column or row vector to be differentiated
%
%    npoints - number of points in the finite difference
%              stencil
%
%    order   - order of the derivative required
%
% Outputs:
%
%    dx      - column or row vector with the derivative
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fdvec.m>

function dx=fdvec(x,npoints,order)

% Check consistency
grumble(x,npoints,order)

% Preallocate the answer
dx=zeros(size(x)); x=x(:);

% Compute edges with sided schemes
for n=1:(npoints-1)/2
    w=fdweights(n,1:npoints,order);
    dx(n)=w(end,:)*x(1:npoints);
    dx(end-n+1)=-w(end,end:-1:1)*x((end-npoints+1):end);
end

% Fill in the middle with centered schemes
stencil=((1-npoints)/2):((npoints-1)/2);
w=fdweights(0,stencil,order);
for n=((npoints-1)/2+1):(numel(x)-(npoints-1)/2)
    dx(n)=w(end,:)*x(stencil+n);
end

end

% Consistency enforcement
function grumble(x,npoints,order)
if (~isnumeric(x))||(~isvector(x))
    error('x argument must be a vector.');
end
if (npoints<1)||(order<1)||(mod(npoints,1)~=0)||(mod(order,1)~=0)
    error('stencil size and derivative order must be positive integers.');
end
if numel(x)<3
    error('x must have more than three elements.');
end
if mod(npoints,2)~=1
    error('the number of stencil points must be odd.');
end
if order>=npoints
    error('derivative order must be smaller than the stencil size.');
end 
end

% So whereas back then I wrote about the tyranny of the majority, today
% I'd combine that with the tyranny of the minorities. [...] I say to
% both bunches, whether you're a majority or minority, bug off! To hell
% with anybody who wants to tell me what to write. All this political
% correctness that's rampant on campuses is bullshit.
%
% You can't fool around with the dangerous notion of telling a university
% what to teach and what not to. If you don't like the curriculum, go to
% another school. Faculty members who toe the same line are sanctimonious
% nincompoops! [...] In the same vein, we should immediately bar all quo-
% tas, which politicize the process through lowered admission standards
% that accept less-qualified students. The terrible result is the price-
% less chance lost by all. The whole concept of higher education is nega-
% ted unless the sole criterion used to determine if students qualify is
% the grades they score on standardized tests.
%
% Ray Bradbury

