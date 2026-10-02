% The function computes the m'th derivative
% Fourier spectral differentiation matrix on grid with N equispa-
% ced points in [0,2pi). Syntax:
%
%                        [x,DM]=fourdif(spin_system,N,m)
% 
% Parameters:
%
%     spin_system - Spinach system with the sys.enable option list
%
%     N - dimension of differentiation matrix
%
%     m - derivative order
%
%  Outputs:
%
%     x - equispaced points 0, 2*pi/N, 4*pi/N, ... , (N-1)*2*pi/N
%
%    DM - m-th order differentiation matrix, or an implicit polyadic
%         when spin_system.sys.enable contains the polyadic option
%
%  Explicit formulae are used to compute the matrices for m=1 and
%  m=2. A discrete Fourier approach is employed for m>2. The prog-
%  ram computes the first column and first row and then uses the 
%  toeplitz() function to create the matrix.
%
%  For m=1 and 2 the code implements a "flipping trick" to improve
%  accuracy as suggested in http://dx.doi.org/10.1137/0916073 
%
%  The polyadic option uses inverse FFT, a numeric CPU diagonal
%  multiplier, and FFT factors with normalised adjoints. On even
%  grids, odd derivatives annihilate the Nyquist mode; even ones
%  retain it. Scale DM by (2*pi/extent)^m for another period. Upload
%  the polyadic once with gpuArray for GPU actions; do not inflate it.
%
%  S.C. Reddy
%  J.A.C. Weideman
%
% <https://spindynamics.org/wiki/index.php?title=fourdif.m>

function [x,DM]=fourdif(spin_system,N,m)

% Check consistency
grumble(spin_system,N,m);

% Run the blob
x=2*pi*(0:N-1)'/N;                           % grid points

% Differentiation matrix
if (nargout>1)&&ismember('polyadic',spin_system.sys.enable)

    % Build the spectral multiplier once on CPU
    mult=fftdiff(m,N,2*pi/N).';
    if (mod(N,2)==0)&&(mod(m,2)==1)
        mult(N/2+1)=0;
    end
    mult=spdiags(mult,0,N,N);

    % Define transforms and their normalised adjoints
    fft_core=struct('action',@(x)fft(x,[],1),...
                    'adjoint',@(x)N*ifft(x,[],1),'dims',[N N]);
    ifft_core=struct('action',@(x)ifft(x,[],1),...
                     'adjoint',@(x)fft(x,[],1)/N,'dims',[N N]);

    % Compose inverse FFT, multiplication, and FFT factors
    DM=polyadic({{ifft_core}})*polyadic({{mult}})*polyadic({{fft_core}});

elseif nargout>1
    h=2*pi/N;                                % grid spacing
    kk=(1:N-1)';
    n1=floor((N-1)/2); 
    n2=ceil((N-1)/2);
    if m==0                                  % compute first column
        col1=[1; zeros(N-1,1)];              % of zeroth derivative
        row1=col1;                           % matrix, which is identity
    elseif m==1                              % compute first column
        if rem(N,2)==0                       % of 1st derivative matrix
            topc=cot((1:n2)'*h/2);    
            col1=[0; 0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]];  
        else
            topc=csc((1:n2)'*h/2);
            col1=[0; 0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
        end
        row1=-col1;                          % first row
    elseif m==2                              % compute first column
        if rem(N,2)==0                       % of 2nd derivative matrix
            topc=csc((1:n2)'*h/2).^2;
            col1=[-pi^2/3/h^2-1/6; -0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
        else           
            topc=csc((1:n2)'*h/2).*cot((1:n2)'*h/2);     
            col1=[-pi^2/3/h^2+1/12; -0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]];
        end        
        row1=col1;                           % first row
    else                                     % employ FFT to compute
        N1=floor((N-1)/2);                   % 1st column of matrix for m>2
        N2=(-N/2)*rem(m+1,2)*ones(rem(N+1,2));
        mwave=1i*[(0:N1) N2 (-N1:-1)];
        col1=real(ifft((mwave.^m).*fft([1 zeros(1,N-1)])));
        if rem(m,2)==0
            row1=col1;                       % first row even derivative
        else
            col1=[0 col1(2:N)]';
            row1=-col1;                      % first row odd derivative
        end
    end   
    DM=toeplitz(col1,row1);
end 
end

% Consistency enforcement
function grumble(spin_system,N,m)
if (~isstruct(spin_system))||(~isfield(spin_system,'sys'))||...
   (~isfield(spin_system.sys,'enable'))||(~iscell(spin_system.sys.enable))
    error('spin_system.sys.enable must be a cell array of option names.');
end
if (~isnumeric(N))||(~isreal(N))||(numel(N)~=1)||(~isfinite(N))||(N<1)||(mod(N,1)~=0)
    error('N must be a positive real integer');
end
if (~isnumeric(m))||(~isreal(m))||(numel(m)~=1)||(~isfinite(m))||(m<1)||(mod(m,1)~=0)
    error('m must be a positive real integer');
end
end

% When you first kissed me, the desire to kiss you back 
% has only very narrowly overpowered the desire to slap 
% you on the face.
%
% One of IK's girlfriends


