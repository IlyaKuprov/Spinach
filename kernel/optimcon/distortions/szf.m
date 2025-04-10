% Applies a discrete single-zero filter:
%
%              Y(k)=X(k)/(1-z)-z*X(k-1)/(1-z);
%
% to a Spinach optimal control module waveform. Treats odd 
% rows of multi-row waveform arrays as real, and even rows
% as imaginary, components of a complex signal. Syntax:
%       
%                      [w,J]=szf(w,z)
%
% Parameters:
%
%    w   - waveform, one time slice per column, and
%          rows arranged as XYXY... with respect to
%          in-phase and quadrature parts on each 
%          control channel
%
%    z   - a vector (one element per XY control pair)
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
% <https://spindynamics.org/wiki/index.php?title=szf.m>

function [w,J]=szf(w,z)

% Check consistency
grumble(w,z);

% Autodiff wrapper
if nargout<2
    
    % Plain call
    w=distort(w(:),z,size(w));

else

    % Autodiff call including Jacobian
    [w,J]=dlfeval(@distort,dlarray(w(:)),z,size(w));

    % Strip autodiff rigging; kill Wirtinger terms
    w=extractdata(w); J=extractdata(J); J=real(J);

end

end

% Actual distortion function
function [w_dist,J]=distort(w,z,dims)

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
    [I,J]=ndgrid(1:ncols,1:ncols);    % Index arrays
    A=(I==1)&(J==1);                  % First element unchanged
    B=((I>1)&(I==J)).*(1/(1-z(n)));   % Main diagonal subsequently
    C=((I>1)&(I==(J+1))).*(-z/(1-z)); % Subdiagonal subsequently

    % Apply the filter
    y=(A+B+C)*transpose(x);

    % Assign back to w_dist
    w_dist(2*n-1,:)=real(y);
    w_dist(2*n,:)=imag(y);

end

% Compute the autodiff Jacobian
if nargout>1, J=dljacobian(w_dist(:),w,1); end

end

% Consistency enforcement
function grumble(w,z)
if (~isnumeric(w))||(~isreal(w))
    error('w must be an array of real numbers.');
end
if mod(size(w,1),2)~=0
    error('the number of rows in w must be even');
end
if (~isnumeric(z))||(~isvector(z))
    error('z must be a vector or a scalar.');
end
if numel(z)~=size(w,1)/2
    error('z must have one element per X,Y control pair.');
end
end

% Bessie Braddock MP: "Winston, you are drunk, and what's more 
%                      you are disgustingly drunk."
%
% Winston Churchill:  "Bessie, my dear, you are ugly, and what's
%                      more, you are disgustingly ugly. But to-
%                      morrow I shall be sober and you will 
%                      still be disgustingly ugly."
%
% Personal communication of Ronald Golding,
% Churchill's bodyguard, related by Richard
% M. Langworth

