% Applies an FIR convolution filter to a Spinach optimal control
% module waveform. Treats odd rows of multi-row waveform arrays
% as real, and even rows as imaginary, components of a complex
% signal. The distal end of the convolution is truncated so the
% output has the same number of samples as the input. Syntax:
%       
%                      [w,J]=firf(w,ker)
%
% Parameters:
%
%    w     - waveform, one time slice per column, and
%            rows arranged as XYXY... with respect to
%            in-phase and quadrature parts on each
%            control channel
%
%    ker   - a vector of FIR filter coefficients
%
% Outputs:
%
%    w     - distorted waveform, same dimension as the
%            input waveform; leaving sufficient ring-
%            down margin is the user's responsibility
%
%    J     - Jacobian matrix with respect to vectorisa-
%            tions of the output and the input arrays
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=firf.m>

function [w,J]=firf(w,ker)

% Check consistency
grumble(w,ker);

% Autodiff wrapper
if nargout<2
    
    % Plain call
    w=distort(w(:),ker,size(w));

else

    % Autodiff call including Jacobian
    [w,J]=dlfeval(@distort,dlarray(w(:)),ker,size(w));

    % Strip autodiff rigging; kill Wirtinger terms
    w=extractdata(w); J=extractdata(J); J=real(J);

end

end

% Actual distortion function
function [w_dist,J]=distort(w,ker,dims)

% Fold into physical dimensions
inp=reshape(w,dims); nrows=dims(1);
ncols=dims(2); nchannels=nrows/2;

% Preallocate the output
w_dist=zeros(size(inp),'like',inp);

% Force coefficients into a column
ker=ker(:);

% Loop over channels
for n=1:nchannels

    % Build complex input signal
    x=inp(2*n-1,:)+1i*inp(2*n,:);

    % Build the convolution matrix
    if ncols>=numel(ker)
        ker_col=[ker; zeros(ncols-numel(ker),1)];
    else
        ker_col=ker(1:ncols);
    end
    ker_row=[ker(1) zeros(1,ncols-1)];
    conv_mat=toeplitz(ker_col,ker_row);

    % Apply the filter
    y=conv_mat*transpose(x);

    % Assign back to w_dist
    w_dist(2*n-1,:)=real(y);
    w_dist(2*n,:)=imag(y);

end

% Compute the autodiff Jacobian
if nargout>1, J=dljacobian(w_dist(:),w,1); end

end

% Consistency enforcement
function grumble(w,ker)
if (~isnumeric(w))||(~isreal(w))
    error('w must be an array of real numbers.');
end
if mod(size(w,1),2)~=0
    error('the number of rows in w must be even.');
end
if (~isnumeric(ker))||(~isvector(ker))
    error('ker must be a vector.');
end
if isempty(ker)
    error('ker must have at least one element.');
end
end

% In theory there is no difference between theory and
% practice. In practice there is.
%
% Yogi Berra
