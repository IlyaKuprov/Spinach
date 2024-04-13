% Checks if a Fokker-Planck state vector has any spatial frequencies
% that its spatial grid is dangerously close to misrepresenting due
% to insufficient point count. Syntax:
%
%                   overwound(rho,spc_dim,spn_dim)
%
% Parameters:
%
%     rho      - Fokker-Planck state vector or
%                a bookshelf stack thereof
%
%     spc_dim  - spatial dimensions of the Fokker-
%                Planck problem, [X Y Z]
%
%     spn_dim  - spin dimension of the Fokker-
%                Planck problem
%
% Output:
%
%     figures and diagnostic messages to the console
%
% Note: if you are running spatial dynamics, such as diffusion and
%       and flow, with finite difference derivative operators, the
%       spatial grid point count shuld be set to several times the
%       minimum Nyquist value.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=overwound.m>

function overwound(rho,spc_dim,spn_dim)

% Check consistency
grumble(rho,spc_dim,spn_dim);

% Find out the stack size
stack_size=size(rho,2);

% Fold the state vector
rho=reshape(rho,[spn_dim spc_dim stack_size]);

% Get spin dimensions out of the way
rho=permute(rho,[2 3 4 1 5]);

% Check X dimension
if spc_dim(1)>1
    
    % Few points are not good
    if spc_dim(1)<10
        error('you must increase point count in X dimension.');
    end
    
    % Run Fourier transform
    rho_kx=fftshift(fft(rho,[],1),1);
    
    % Pull X dimension forward
    rho_kx=permute(rho_kx,[1 2 3 4 5]);
    rho_kx=reshape(rho_kx,[size(rho_kx,1) ...
                           numel(rho_kx)/size(rho_kx,1)]);
    
    % Draw the plot
    freq_axis=linspace(-1,1,size(rho_kx,1));
    figure(); plot(freq_axis,sum(abs(rho_kx),2));
    xlabel('X spatial frequency relative to Nyquist limit');
    ylabel('population density, a.u.'); kgrid;
    
end

% Check Y dimension
if spc_dim(2)>1
    
    % Few points are not good
    if spc_dim(2)<10
        error('you must increase point count in Y dimension.');
    end
    
    % Run Fourier transform
    rho_ky=fftshift(fft(rho,[],2),2);
    
    % Pull Y dimension forward
    rho_ky=permute(rho_ky,[2 1 3 4 5]);
    rho_ky=reshape(rho_ky,[size(rho_ky,1) ...
                           numel(rho_ky)/size(rho_ky,1)]);
    
    % Draw the plot
    freq_axis=linspace(-1,1,size(rho_ky,1));
    figure(); plot(freq_axis,sum(abs(rho_ky),2));
    xlabel('Y spatial frequency relative to Nyquist limit');
    ylabel('population density, a.u.'); kgrid;
    
end

% Check X dimension
if spc_dim(3)>1
    
    % Few points are not good
    if spc_dim(3)<10
        error('you must increase point count in Z dimension.');
    end
    
    % Run Fourier transform
    rho_kz=fftshift(fft(rho,[],3),3);
    
    % Pull X dimension forward
    rho_kz=permute(rho_kz,[3 2 1 4 5]);
    rho_kz=reshape(rho_kz,[size(rho_kz,1) ...
                           numel(rho_kz)/size(rho_kz,1)]);
    
    % Draw the plot
    freq_axis=linspace(-1,1,size(rho_kz,1));
    figure(); plot(freq_axis,sum(abs(rho_kz),2));
    xlabel('Z spatial frequency relative to Nyquist limit');
    ylabel('population density, a.u.'); kgrid;
 
end

end

% Consistency enforcement
function grumble(rho,spc_dim,spn_dim)
if (~isnumeric(rho))
    error('rho must be a numeric array.');
end
if (~isnumeric(spc_dim))||(~isreal(spc_dim))||(~isvector(spc_dim))||...
   (any(spc_dim<1))||(any(mod(spc_dim,1)~=0))||(numel(spc_dim)~=3)
    error('spc_dim must be a vector of three positive integers.');
end
if (~isnumeric(spn_dim))||(~isreal(spn_dim))||(~isscalar(spn_dim))||...
   (spn_dim<1)||(mod(spn_dim,1)~=0)
    error('spn_dim must be a positive integer.');
end
if size(rho,1)~=prod([spn_dim spc_dim])
    error('the number of rows in rho is not consistent with dimensions specified.');
end
end

% "Why quantum mechanics? And why turbulence?" - the two questions 
% Werner Heisenberg wanted to ask God. He was quite sure God could
% answer the first one.

