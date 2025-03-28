% Applies a discrete single-pole filter:
%
%                 Y(n)=(1-p)*X(n)+p*Y(n-1)
%
% to a Spinach optimal control module waveform. Treats odd 
% rows of multi-row waveform arrays as real, and even rows
% as imaginary, components of a complex signal. Syntax:
%       
%                      [w,J]=spf(w,p)
%
% Parameters:
%
%    w   - waveform, one time slice per column, and
%          rows arranged as XYXY... with respect to
%          in-phase and quadrature parts on each 
%          control channel
%
%    p   - a vector (one element per XY control pair)
%          containing the filter coefficient:
%
%             p=exp(-r*dt+1i*(omega-omega_rf)*dt)
%
%          where r is th damping rate, omega is the
%          pole frequency, omega_rf is the rotating
%          frame frequency, and dt is the time dis-
%          cretisation step.
%
% Outputs:
%
%    w   - distorted waveform, same dimension as the
%          input waveform; leaving sufficient ring-
%          down margin is the user's responsibility
%
%    J   - Jacobian matrix with respect to vectorisa-
%          tions of the output and the input arrays
%
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=spf.m>

function [w,J]=spf(w,p)

% Check consistency
grumble(w,p);

% Autodiff wrapper
if nargout<2
    
    % Plain call
    w=distort(w(:),p,size(w));

else

    % Autodiff call including Jacobian
    [w,J]=dlfeval(@distort,dlarray(w(:)),p,size(w));

    % Strip autodiff rigging; kill Wirtinger terms
    w=extractdata(w); J=extractdata(J); J=real(J);

end

end

% Actual distortion function
function [w_dist,J]=distort(w,p,dims)

% Fold into physical dimensions
inp=reshape(w,dims); nrows=dims(1); 
ncols=dims(2); nchannels=nrows/2;

% Preallocate the output
w_dist=zeros(size(inp),'like',inp);

% Loop over channels
for n=1:nchannels

    % Build complex input signal
    x=inp(2*n-1,:)+1i*inp(2*n,:);

    % Build the filter matrix
    [I,J]=ndgrid(1:ncols,1:ncols);             % Index arrays
    A=(J==1).*p(n).^(I-1);                     % First element unchanged
    B=((J>=2)&(I>=J)).*(1-p(n)).*p(n).^(I-J);  % The rest is filtered

    % Apply the filter
    y=(A+B)*transpose(x);

    % Assign back to w_dist
    w_dist(2*n-1,:)=real(y);
    w_dist(2*n,:)=imag(y);

end

% Compute the autodiff Jacobian
if nargout>1, J=dljacobian(w_dist(:),w,1); end

end

% Consistency enforcement
function grumble(w,p)
if (~isnumeric(w))||(~isreal(w))
    error('w must be an array of real numbers.');
end
if mod(size(w,1),2)~=0
    error('the number of rows in w must be even.');
end
if (~isnumeric(p))||(~isvector(p))
    error('p must be a vector or a scalar.');
end
if numel(p)~=size(w,1)/2
    error('p must have one element per X,Y control pair.');
end
if any(abs(p(:))>=1)
    error('must have |p(k)|<1 for causal stability.');
end
end

% All people serve their ambition. In that matter, there 
% are no atheists. There are only people who know, and 
% don't know, what God they serve.
%
% Jordan Peterson, 
% 12 Rules for Life: an Antidote to Chaos

