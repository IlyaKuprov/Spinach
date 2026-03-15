% DPFGSE signal selection, based on Equation 3 from the paper by
% Stott et al. (https://doi.org/10.1006/jmre.1997.1110). Syntax:
%
%      fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)
%
% Parameters:
%
%    parameters.g_amp       - amplitudes of the two gradients, T/m
%
%    parameters.g_dur       - gradient duration, seconds
%
%    parameters.rf_frq_list - soft pulse parameters that will
%    parameters.rf_amp_list   be passed to shaped_pulse_af
%    parameters.rf_dur_list   function
%    parameters.rf_phi
%    parameters.max_rank
%
%    parameters.sweep       - detection sweep width, Hz
%
%    parameters.npoints     - number of points in the fid
%
%    parameters.offset      - transmitter and receiver offset, Hz
%
% Outputs:
%
%    fid - free induction decay of what is effectively a
%          1D pulse-acquire NMR experiment
%
% Notes: at least a hundred points are required in the spatial
%        dimension; increase until the answer stops changing.
%
% i.kuprov@soton.amann.ac.uk
% p.lally@soton.ac.uk
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dpfgse_select.m>

function fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(parameters);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+','1H');
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Hard 90 on everything
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Gradient
rho=step(spin_system,L+parameters.g_amp(1)*G{1},rho,parameters.g_dur);
 
% Soft 180 on user-specified frequency
rho=shaped_pulse_af(spin_system,L,Lx,Ly,rho,parameters.rf_frq_list-parameters.offset,...
                    parameters.rf_amp_list,parameters.rf_dur_list,...
                    parameters.rf_phi,parameters.max_rank,'expv');
 
% Gradients
rho=step(spin_system,L+parameters.g_amp(1)*G{1},rho,parameters.g_dur);
rho=step(spin_system,L+parameters.g_amp(2)*G{1},rho,parameters.g_dur);
 
% Soft 180 on user-specified frequency
rho=shaped_pulse_af(spin_system,L,Lx,Ly,rho,parameters.rf_frq_list-parameters.offset,...
                    parameters.rf_amp_list,parameters.rf_dur_list,...
                    parameters.rf_phi,parameters.max_rank,'expv');
  
% Gradient
rho=step(spin_system,L+parameters.g_amp(2)*G{1},rho,parameters.g_dur);

% Run the evolution and watch the coil state
fid=evolution(spin_system,L,parameters.coil,rho,...
    1/parameters.sweep,parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'rho0')
    error('parameters.rho0 field must be present.');
end
if ~isnumeric(parameters.rho0)
    error('parameters.rho0 must be numeric.');
end
if ~isfield(parameters,'coil')
    error('parameters.coil field must be present.');
end
if ~isnumeric(parameters.coil)
    error('parameters.coil must be numeric.');
end
if ~isfield(parameters,'npts')
    error('parameters.npts field must be present.');
end
if (~isnumeric(parameters.npts))||(~isreal(parameters.npts))||...
   (~isvector(parameters.npts))||any(parameters.npts<1)||...
   any(mod(parameters.npts,1)~=0)
    error('parameters.npts must be a vector of positive integers.');
end
if ~isfield(parameters,'offset')
    error('parameters.offset field must be present.');
end
if (~isnumeric(parameters.offset))||(~isreal(parameters.offset))||(~isscalar(parameters.offset))
    error('parameters.offset must be a real scalar.');
end
if ~isfield(parameters,'rf_frq_list')
    error('parameters.rf_frq_list field must be present.');
end
if (~isnumeric(parameters.rf_frq_list))||(~isreal(parameters.rf_frq_list))||...
   (~isvector(parameters.rf_frq_list))
    error('parameters.rf_frq_list must be a real vector.');
end
if ~isfield(parameters,'rf_amp_list')
    error('parameters.rf_amp_list field must be present.');
end
if (~isnumeric(parameters.rf_amp_list))||(~isreal(parameters.rf_amp_list))||...
   (~isvector(parameters.rf_amp_list))
    error('parameters.rf_amp_list must be a real vector.');
end
if ~isfield(parameters,'rf_dur_list')
    error('parameters.rf_dur_list field must be present.');
end
if (~isnumeric(parameters.rf_dur_list))||(~isreal(parameters.rf_dur_list))||...
   (~isvector(parameters.rf_dur_list))||any(parameters.rf_dur_list<=0)
    error('parameters.rf_dur_list must be a positive real vector.');
end
if ~isfield(parameters,'rf_phi')
    error('parameters.rf_phi field must be present.');
end
if (~isnumeric(parameters.rf_phi))||(~isreal(parameters.rf_phi))||(~isscalar(parameters.rf_phi))
    error('parameters.rf_phi must be a real scalar.');
end
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank field must be present.');
end
if (~isnumeric(parameters.max_rank))||(~isreal(parameters.max_rank))||...
   (~isscalar(parameters.max_rank))||(parameters.max_rank<1)||...
   (mod(parameters.max_rank,1)~=0)
    error('parameters.max_rank must be a positive integer.');
end
if ~isfield(parameters,'sweep')
    error('parameters.sweep field must be present.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isscalar(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('parameters.npoints field must be present.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (~isscalar(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a positive integer.');
end
if ~isfield(parameters,'g_amp')
    error('gradient amplitudes should be provided in parameters.g_amp field.');
end
if (~isnumeric(parameters.g_amp))||(~isreal(parameters.g_amp))||...
   (~isvector(parameters.g_amp))||(numel(parameters.g_amp)~=2)
    error('parameters.g_amp must be a vector with two real numbers.');
end
if ~isfield(parameters,'g_dur')
    error('gradient duration should be provided in parameters.g_dur field.');
end
if (~isnumeric(parameters.g_dur))||(~isreal(parameters.g_dur))||...
   (~isscalar(parameters.g_dur))||(parameters.g_dur<=0)
    error('parameters.g_dur must be a positive real number.');
end
if (numel(parameters.rf_frq_list)~=numel(parameters.rf_amp_list))||...
   (numel(parameters.rf_amp_list)~=numel(parameters.rf_dur_list))
    error('parameters.rf_frq_list, parameters.rf_amp_list, and parameters.rf_dur_list must have the same number of elements.');
end
end

% We (Mr Rosen and I) had sent you our manuscript for publication and 
% had not authorised you to show it to specialists before it is print-
% ed. I see no reason to address the - in any case erroneous - comments
% of your anonymous expert. On the basis if this incident I prefer to
% publish the paper elsewhere.
% 
% Albert Einstein, to John Tate (the editor of Physical Review),
% in 1936, about a paper titled "Do Gravitational Waves Exist?"

