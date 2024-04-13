% Cost function for optimal control using the GRAPE algorithm. Returns 
% fidelity, gradient and Hessian for a given waveform, specified in po-
% lar coordinates. Only the phase channel gradient is returned, the am-
% plitude profile is taken as a given. Syntax:
%
%   [traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)
%
% Parameters:
%
%   phi_profile   - set of control pulse phases from an amplitude-phase 
%                   description.
%
% Outputs:
%
%   fidelity      - figure of merit for the overlap of the current state
%                   of the system and the desired state(s). When penalty
%                   methods are specified, fidelity is returned as an ar-
%                   ray separating the penalties from the simulation
%                   fidelity.
%
%   gradient      - gradient of the fidelity with respect to the control 
%                   sequence. When penalty methods are specified, gradi-
%                   ent is returned as an array separating penalty gra-
%                   dients from the fidelity gradient.
%
%   hessian       - Hessian of the fidelity with respect to the control 
%                   sequence. When penalty methods are specified, gradi-
%                   ent is returned as an array separating penalty Hes-
%                   sians from the fidelity Hessian.
%
% david.goodwin@inano.au.dk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grape_phase.m>

function [traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)

% Check consistency
grumble(spin_system,spin_system.control.amplitudes,phi_profile);

% Count penalty terms
npenterms=numel(spin_system.control.penalties);

% Transform into Cartesian coordinates
waveform_xy=zeros(2*size(phi_profile,1),size(phi_profile,2));
for n=1:size(phi_profile,1)
    [waveform_xy(2*n-1,:),waveform_xy(2*n,:)]=polar2cartesian(spin_system.control.amplitudes(n,:),phi_profile(n,:));
end

% Just fidelity
if nargout==2

    % Call Cartesian GRAPE
    [traj_data,fidelity]=grape_xy(waveform_xy,spin_system);
    
% Fidelity and gradient
elseif nargout==3
    
    % Call Cartesian GRAPE
    [traj_data,fidelity,grad_xy]=grape_xy(waveform_xy,spin_system);
    
    % Preallocate the answer
    gradient=zeros(size(phi_profile,1),size(phi_profile,2),npenterms+1);
    
    % Loop over phase tracks
    for n=1:size(phi_profile,1)
        
        % Loop over penalites
        for k=1:(npenterms+1)
            
            % Translate derivatives
            [~,~,~,gradient(n,:,k)]=cartesian2polar(waveform_xy(2*n-1,:),waveform_xy(2*n,:),...
                                                    grad_xy(2*n-1,:,k),  grad_xy(2*n,:,k));
                                                
        end
        
    end

% Fidelity, gradient and Hessian
elseif nargout==4
    
    % Call Cartesian ensemble GRAPE
    [traj_data,fidelity,grad_xy,hess_xy]=grape_xy(waveform_xy,spin_system);
    
    % Preallocate the answer
    gradient=zeros(size(phi_profile,1),size(phi_profile,2),npenterms+1);
    hessian=zeros(numel(phi_profile),numel(phi_profile),npenterms+1);
    
    % Loop over phase tracks
    for n=1:size(phi_profile,1)
        
        % Loop over penalites
        for k=1:(npenterms+1)
            
            % Translate derivatives
            [~,~,~,gradient(n,:,k)]=cartesian2polar(waveform_xy(2*n-1,:),waveform_xy(2*n,:),...
                                                    grad_xy(2*n-1,:,k),  grad_xy(2*n,:,k));
                                                
        end
        
    end
    
    % Transform Hessian from n^2 kxk block matrices to k^2 nxn matrices
    for n=1:size(hess_xy,3)
        hess_xy(:,:,n)=hess_reorder(hess_xy(:,:,n),numel(spin_system.control.operators),...
                                                   spin_system.control.pulse_nsteps);
    end
    
    % Loop over phase tracks
    for n=1:size(phi_profile,1)
        for k=1:size(phi_profile,1)
            
            % Waveforms
            fx=waveform_xy(2*n-1,:);
            fy=waveform_xy(2*n,:);
            
            % Loop over penalties
            for m=1:(npenterms+1)
            
                % Gradients
                Dx=grad_xy(2*n-1,:,m);
                Dy=grad_xy(2*n,:,m);
                
                % Hessian blocks
                hess_block=hess_xy(1+2*spin_system.control.pulse_nsteps*(n-1):2*spin_system.control.pulse_nsteps*(n),...
                                   1+2*spin_system.control.pulse_nsteps*(k-1):2*spin_system.control.pulse_nsteps*(k),m);
                Dxx=hess_block(1:end/2,1:end/2);     Dxy=hess_block(1:end/2,1+end/2:end);
                Dyx=hess_block(1+end/2:end,1:end/2); Dyy=hess_block(1+end/2:end,1+end/2:end);
            
                % Translate second derivatives
                [~,~,~,~,~,~,~,Dpp]=cartesian2polar(fx,fy,Dx,Dy,Dxx,Dxy,Dyx,Dyy);
            
                % Make the phase Hessian
                hessian(1+spin_system.control.pulse_nsteps*(n-1):spin_system.control.pulse_nsteps*(n),...
                        1+spin_system.control.pulse_nsteps*(k-1):spin_system.control.pulse_nsteps*(k),m)=Dpp;
                    
            end
                    
        end
        
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,amp_profile,phi_profile)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('optimal control module requires Lioville space formalism.');
end
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if mod(numel(spin_system.control.operators),2)~=0
    error('phase-amplitude optimisations must have an even number of controls.');
end
if (~isnumeric(amp_profile))||(~isreal(amp_profile))
    error('amp_profile must be an array of real numbers.');
end
if any(amp_profile(:)<0)
    error('all elements of amp_profile must be non-negative.');
end
if (~isnumeric(phi_profile))||(~isreal(phi_profile))
    error('phi_profile must be an array of real numbers');
end
if numel(amp_profile)~=numel(phi_profile)
    error('amp_profile and phi_profile must have the same number of elements.');
end
if size(amp_profile,1)~=numel(spin_system.control.operators)/2
    error('the number of rows in amp_profile must be half the number of controls.');
end
if size(phi_profile,1)~=numel(spin_system.control.operators)/2
    error('the number of rows in phi_profile must be half the number of controls.');
end
switch spin_system.control.integrator
    case 'rectangle'
        if size(amp_profile,2)~=spin_system.control.pulse_nsteps
            error('the number of columns in amp_profile must be equal to the number of time steps.');
        end
        if size(phi_profile,2)~=spin_system.control.pulse_nsteps
            error('the number of columns in phi_profile must be equal to the number of time steps.');
        end
    case 'trapezium'
        if size(amp_profile,2)~=(spin_system.control.pulse_nsteps+1)
            error('the number of columns in amp_profile must be (number of time steps)+1.');
        end
        if size(phi_profile,2)~=(spin_system.control.pulse_nsteps+1)
            error('the number of columns in phi_profile must be (number of time steps)+1.');
        end
    otherwise
        error('unknown time propagation algorithm.');
end
end

% Perfection (in design) is achieved not when there is nothing 
% more to add, but rather when there is nothing more to take away.
%
% Antoine de Saint-Exupery

