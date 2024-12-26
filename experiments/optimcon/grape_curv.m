% Cost function for optimal control using the GRAPE algorithm. Returns 
% fidelity and gradient for a given waveform, specified in arbitrary
% curvilinear coordinates. Syntax:
%
%      [traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,...
%                                            dx_du,spin_system)
%
% Parameters:
%
%   waveform_u   -  pulse waveform in curvilinear coordinates with indi-
%                   vidual coordinates in columns and time in rows
%
%   u2x          -  a handle to a function that takes a column of curvi-
%                   linear coordinates and returns a column of coeffici-
%                   ents in front of the control operators
%
%   dx_du        -  a handle to a function that takes a column of curvi-
%                   linear coordinates and returns the Jacobian matrix
%                   with the following structure:
%                    
%                       [dx(1)_du(1) dx(2)_du(1) dx(3)_du(1) ...
%                        dx(1)_du(2) dx(2)_du(2) dx(3)_du(2) ...
%                           ...         ...         ...      ...]
%
% Outputs:
%
%   traj_data     - system trajectory data structure used for visualisa-
%                   tion and progress reports
%
%   fidelity      - figure of merit for the overlap of the current state
%                   of the system and the desired state(s). When penalty
%                   methods are specified, fidelity is returned as an ar-
%                   ray separating the penalties from the simulation
%                   fidelity.
%
%   df_du         - gradient of the fidelity with respect to the control 
%                   sequence. When penalty methods are specified, gradi-
%                   ent is returned as an array separating penalty gra-
%                   dients from the fidelity gradient.
%
% Note: penalities are computed using the rectilinear representation.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=grape_curv.m>

function [traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,...
                                               dx_du,spin_system)
% Check consistency
grumble(waveform_u,u2x,dx_du,spin_system);

% Count rectilinear and curvilinear coordinates
size_x=numel(spin_system.control.operators); 
n_time_pts=size(waveform_u,2);
size_u=size(waveform_u,1); 

% Transform into rectilinear coordinates
waveform_x=zeros(size_x,n_time_pts);
for n=1:n_time_pts
    waveform_x(:,n)=u2x(waveform_u(:,n));
end

% Decide the output
switch nargout
    
    % Just fidelity
    case 2

        % Call rectilinear coordinate GRAPE
        [traj_data,fidelity]=grape_xy(waveform_x,spin_system);
    
    % Fidelity and gradient
    case 3
    
        % Call rectilinear coordinate GRAPE
        [traj_data,fidelity,df_dx]=grape_xy(waveform_x,spin_system);
        
        % Preallocate curvilinear gradients
        df_du=zeros(size_u,n_time_pts,size(df_dx,3));
        
        % Fidelity and penalties
        for n=1:size(df_dx,3)
            
            % Loop over time points
            for k=1:n_time_pts
                
                % Translate gradients into curvilinear coordinates
                df_du(:,k,n)=dx_du(waveform_u(:,k))*df_dx(:,k,n);
                
            end
            
        end
       
    otherwise
        
        % Complain and bomb out
        error('incorrect number of output arguments');
        
end

end

% Consistency enforcement
function grumble(waveform_u,u2x,dx_du,spin_system)
if ~isfield(spin_system,'control')
    error('spin_system object lacks control information, run optimcon() first.');
end
if (~isa(u2x,'function_handle'))||(~isa(dx_du,'function_handle'))
    error('u2x and dx_du must be function handles.');
end
if (~isnumeric(waveform_u))||(~isreal(waveform_u))
    error('waveform_u must be an array of real numbers.');
end
end

% The worst case of unintended consequences IK has ever had was a paper, submitted
% to Science Advances, where a new MRI contrast agent was proposed and MRI images
% reported for animals, including a beagle dog. An important question was about how
% complete was the clearance after a few days - the substance did turn up in the
% urine, but the extent of its retention in the body had not been quantified. With
% his editor hat on, IK took the same position as the Reviewers and made the accep-
% tance conditional on the authors providing some evidence that the contrast agent
% did not accumulate in the organism. Two weeks later, a revised version of the pa-
% per turned up... the authors killed the dog and ran its tissues through ICP-MS. 

