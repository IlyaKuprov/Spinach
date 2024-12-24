% 2D imaging sequence with spiral sampling of the k-space. Syntax:
%
%   mri=spiral_pulse_sequence(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. 
%
% Parameters:
%
%      parameters.t_echo     - echo time, seconds
%
%      parameters.spiral_frq - frequency of the spiral, Hz 
%
%      parameters.spiral_dur - duration of the spiral, seconds
%
%      parameters.grad_amp   - gradient amplitude at the end 
%                              of the spiral, T/m
%
% Outputs:
%
%      mri - MRI image with square sinebell apodisation
%
% Prior to being Fourier transformed, the spiral is resampled onto
% a suitable square grid - change the sequence to output spiral_x,
% spiral_y and spiral_z variables if you plan to process the data
% in a different way.
%
% ilya.kuprov@weizmann.ac.uk
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=spiral.m>

function mri=spiral(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Make pulse operators
Hp=operator(spin_system,'L+','1H');
Hy=kron(speye(prod(parameters.npts)),(Hp-Hp')/2i);

% Apply 90-degree pulse
rho=step(spin_system,Hy,parameters.rho0,pi/2);

% Evolve the for the echo time
rho=evolution(spin_system,L,[],rho,parameters.t_echo,1,'final');

% Apply 180-degree pulse
rho=step(spin_system,Hy,rho,pi);

% Evolve the for the echo time
rho=evolution(spin_system,L,[],rho,parameters.t_echo,1,'final');

% Build the spiral
time_grid=linspace(0,parameters.spiral_dur,parameters.spiral_npts);
spiral_x=time_grid*(parameters.grad_amp/parameters.spiral_dur).*cos(parameters.spiral_frq*time_grid);
spiral_y=time_grid*(parameters.grad_amp/parameters.spiral_dur).*sin(parameters.spiral_frq*time_grid);

% Report the sampling spiral to the user
figure(); plot(spiral_x,spiral_y,'b-'); ktitle('$k$-space sampling spiral');
kxlabel('X gradient amplitude, T/m'); kylabel('Y gradient amplitude, T/m'); 
kgrid; axis equal; drawnow();

% Prelocate the spiral trajectory array
spiral_z=zeros([parameters.spiral_npts 1],'like',1i);

% Get the time step
time_step=(parameters.spiral_dur/parameters.spiral_npts);

% Start feedback timer
feedback=tic();

% Propagate through the spiral
for n=1:parameters.spiral_npts
    
    % Record trajectory point
    spiral_z(n)=parameters.coil'*rho;
    
    % Take a step forward
    rho=step(spin_system,L+spiral_x(n)*G{1}+spiral_y(n)*G{2},rho,time_step);
    
    % Inform the user
    if (n==parameters.spiral_npts)||(toc(feedback)>1)
        report(spin_system,['spiral trajectory step ' num2str(n) '/'...
               num2str(parameters.spiral_npts) '...']);
        feedback=tic();
    end
    
end

% Resample onto square grid
x_grid=linspace(-parameters.grad_amp,parameters.grad_amp,sqrt(parameters.spiral_npts)/2);
y_grid=linspace(-parameters.grad_amp,parameters.grad_amp,sqrt(parameters.spiral_npts)/2);
[X,Y]=ndgrid(x_grid,y_grid); fid=griddata(spiral_x,spiral_y,spiral_z,X,Y,'cubic');

% Replace NaN values with zeros
fid(isnan(fid))=0;

% Apodisation
fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}});

% Fourier transform
mri=real(fftshift(fft2(ifftshift(fid))));

% Transpose and flip around
mri=fliplr(mri');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,G,F)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function is only available in sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~isnumeric(F))||(~ismatrix(H))||(~ismatrix(R))||...
   (~ismatrix(K))||(~ismatrix(F))
    error('H,R,K,F arguments must be matrices.');
end
if (~all(size(H)==size(R)))||...
   (~all(size(R)==size(K)))||...
   (~all(size(K)==size(F)))
    error('H,R,K,F matrices must have the same dimension.');
end
if ~iscell(G)
    error('the G argument must be a cell array.');
end
if ~isfield(parameters,'t_echo')
    error('echo time must be specified in parameters.t_echo field.');
end
if (~isnumeric(parameters.t_echo))||(~isreal(parameters.t_echo))||...
   (~isscalar(parameters.t_echo))||(parameters.t_echo<=0)
    error('parameters.t_echo must be a positive real scalar.');
end
if ~isfield(parameters,'spiral_frq')
    error('spiral frequency must be specified in parameters.spiral_frq field.');
end
if (~isnumeric(parameters.spiral_frq))||(~isreal(parameters.spiral_frq))||...
   (~isscalar(parameters.spiral_frq))
    error('parameters.spiral_frq must be a real scalar.');
end
if ~isfield(parameters,'spiral_dur')
    error('spiral duration must be specified in parameters.spiral_dur field.');
end
if (~isnumeric(parameters.spiral_dur))||(~isreal(parameters.spiral_dur))||...
   (~isscalar(parameters.spiral_dur))||(parameters.spiral_dur<=0)
    error('parameters.spiral_dur must be a positive real scalar.');
end
if ~isfield(parameters,'grad_amp')
    error('gradient amplitude must be specified in parameters.grad_amp field.');
end
if (~isnumeric(parameters.grad_amp))||(~isreal(parameters.grad_amp))||...
   (~isscalar(parameters.grad_amp))
    error('parameters.grad_amp must be a real scalar.');
end
end

% According to a local legend at Oxford's Magdalen College, during the
% the Second World War meat shortages, the College's deer were reclass-
% ified as vegetables to avoid requisition by the Ministry of Food.

