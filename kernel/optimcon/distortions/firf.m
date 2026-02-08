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

% Force coefficients into a column
ker=ker(:);

% Fold into physical dimensions
dims=size(w); nrows=dims(1);
ncols=dims(2); nchannels=nrows/2;

% Build the convolution matrix
if ncols>=numel(ker)
    ker_col=[ker; zeros(ncols-numel(ker),1)];
else
    ker_col=ker(1:ncols);
end
ker_row=[ker(1) zeros(1,ncols-1)];
conv_mat=toeplitz(ker_col,ker_row);
conv_mat=sparse(conv_mat);

% Preallocate the output
w_dist=zeros(size(w),'like',w);

% Precompute Jacobian components
if nargout>1

    % Build sparse identity matrix
    id_cols=speye(ncols);

    % Separate real and imaginary parts
    a_mat=real(conv_mat);
    b_mat=imag(conv_mat);

    % Preallocate Jacobian matrix
    J=spalloc(numel(w),numel(w),4*nchannels*nnz(conv_mat));

end

% Loop over channels
for n=1:nchannels

    % Build complex input signal
    x=w(2*n-1,:)+1i*w(2*n,:);

    % Apply the filter
    y=conv_mat*transpose(x);

    % Assign back to w_dist
    w_dist(2*n-1,:)=real(y);
    w_dist(2*n,:)=imag(y);

    % Compute Jacobian block
    if nargout>1

        % Build row selectors
        row_re=2*n-1; row_im=2*n;
        sel_re=kron(id_cols,sparse(1,row_re,1,1,nrows));
        sel_im=kron(id_cols,sparse(1,row_im,1,1,nrows));

        % Build row inserters
        ins_re=kron(id_cols,sparse(row_re,1,1,nrows,1));
        ins_im=kron(id_cols,sparse(row_im,1,1,nrows,1));

        % Accumulate Jacobian contribution
        J=J+ins_re*(a_mat*sel_re-b_mat*sel_im)+...
            ins_im*(b_mat*sel_re+a_mat*sel_im);

    end

end

% Return distorted waveform
w=w_dist;

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

% In 1952, claims that smoking causes cancer led Kent
% Cigarettes to introduce an asbestos filter to "pro-
% tect" smokers. They called it "micronite filter".

