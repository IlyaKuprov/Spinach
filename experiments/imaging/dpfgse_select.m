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
% i.kuprov@soton.ac.uk
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
if ~isfield(parameters,'g_amp')
    error('gradient amplitudes should be provided in parameters.g_amp field.');
end
if (~isreal(parameters.g_amp))||(numel(parameters.g_amp)~=2)
    error('parameters.g_amp must be a vector with two real numbers.');
end
if ~isfield(parameters,'g_dur')
    error('gradient duration should be provided in parameters.g_dur field.');
end
if (~isreal(parameters.g_dur))||(~isscalar(parameters.g_dur))||...
   (parameters.g_dur<=0)
    error('parameters.g_dur must be a positive real number.');
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

