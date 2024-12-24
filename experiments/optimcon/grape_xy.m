% Cost function for optimal control using the GRAPE algorithm. Returns 
% fidelity, gradient and hessian for a given waveform, specified in Car-
% tesian coordinates (x and y channels). Syntax:
%
%    [traj_data,fidelity,grad,hess]=grape_xy(waveform,spin_system)
%
% Parameters:
%
%   waveform      - normalised set of control amplitudes.
%
% Outputs:
%
%   traj_data     - trajectory data
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
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grape_xy.m>

function [traj_data,fidelity,grad,hess]=grape_xy(waveform,spin_system)

% Check consistency
grumble(spin_system,waveform);

% Count penalty terms
npenterms=numel(spin_system.control.penalties);

% Translate the basis if necessary
if ~isempty(spin_system.control.basis)
    nwaves=size(spin_system.control.basis,1);
    waveform=waveform*spin_system.control.basis;
end

% Decide how much needs computing
if nargout==2
    
    % Preallocate output
    fidelity=zeros(1,npenterms+1);
    
    % Calculate the objective and its gradient
    [traj_data,fidelity(1)]=ensemble(waveform,spin_system);
    
    % Apply all penalties
    for n=1:npenterms
        
        % Calculate the penalty (this should be moved to ensemble)
        pen=penalty(waveform,spin_system.control.penalties{n},...
                             spin_system.control.l_bound,...
                             spin_system.control.u_bound);
        
        % Add to fidelity array
        fidelity(n+1)=spin_system.control.p_weights(n)*pen;
        
    end
    
elseif nargout==3
    
    % Preallocate output
    fidelity=zeros(1,npenterms+1);
    grad=zeros(size(waveform,1),size(waveform,2),npenterms+1);
    
    % Calculate the objective and its gradient
    [traj_data,fidelity(1),grad(:,:,1)]=ensemble(waveform,spin_system);
    
    % Apply all penalties
    for n=1:npenterms
        
        % Calculate the penalty
        [pen,pen_grad]=penalty(waveform,spin_system.control.penalties{n},...
                                        spin_system.control.l_bound,...
                                        spin_system.control.u_bound);
        
        % Add to relevant arrays
        fidelity(n+1)=spin_system.control.p_weights(n)*pen; 
        grad(:,:,n+1)=spin_system.control.p_weights(n)*pen_grad;
        
    end
    
    % Translate the basis if necessary
    if ~isempty(spin_system.control.basis)
        grad_in_basis=zeros(size(waveform,1),nwaves,npenterms+1);
        for n=1:(npenterms+1)
            grad_in_basis(:,:,n)=grad(:,:,n)*spin_system.control.basis';
        end
        grad=grad_in_basis;
    end
    
elseif nargout==4
    
    % Preallocate output
    fidelity=zeros(1,npenterms+1);
    grad=zeros(size(waveform,1),size(waveform,2),npenterms+1);
    hess=zeros(numel(waveform),numel(waveform),npenterms+1);
    
    % Calculate objective, gradient and Hessian
    [traj_data,fidelity(1),grad(:,:,1),hess(:,:,1)]=ensemble(waveform,spin_system);
    
    % Apply all penalties
    for n=1:npenterms
        
        % Calculate the penalty
        [pen,pen_grad,pen_hess]=penalty(waveform,spin_system.control.penalties{n},...
                                                 spin_system.control.l_bound,...
                                                 spin_system.control.u_bound);
        
        % Add to relevant arrays
        fidelity(:,n+1)=spin_system.control.p_weights(n)*pen; 
        grad(:,:,n+1)=spin_system.control.p_weights(n)*pen_grad; 
        hess(:,:,n+1)=spin_system.control.p_weights(n)*pen_hess;
        
    end
    
    % Translate the basis if necessary
    if ~isempty(spin_system.control.basis)
        
        % Transform the gradient
        grad_in_basis=zeros(size(waveform,1),nwaves,npenterms+1);
        for n=1:(npenterms+1)
            grad_in_basis(:,:,n)=grad(:,:,n)*spin_system.control.basis';
        end
        grad=grad_in_basis;
        
        % Preallocate the transform
        hess_in_basis=zeros([numel(spin_system.control.operators)*nwaves ...
                             numel(spin_system.control.operators)*nwaves npenterms+1]);
        
        % Transform Hessian and penalties
        for n=1:size(hess,3)
            
            % Reorder Hessian
            hess_re=hess_reorder(hess(:,:,n),numel(spin_system.control.operators),...
                                                spin_system.control.pulse_nsteps);
            
            % Unfold block matrix to 4d tensor
            hess_re=permute(reshape(hess_re,[spin_system.control.pulse_nsteps ...
                                             numel(spin_system.control.operators) ...
                                             spin_system.control.pulse_nsteps ...
                                             numel(spin_system.control.operators)]),[1 3 2 4]);
            % Preallocate the transform
            hess_trans=zeros([nwaves ...
                              nwaves ...
                              numel(spin_system.control.operators) ...
                              numel(spin_system.control.operators)]);
            
            % Compute scalar products
            for j=1:numel(spin_system.sys.controls)
                for k=1:numel(spin_system.sys.controls)
                    hess_trans(:,:,j,k)=spin_system.control.basis*...
                                       (squeeze(hess_re(:,:,j,k))*...
                                        spin_system.control.basis');
                end
            end
            
            % Fold 4d tensor to block matrix
            hess_trans=reshape(permute(hess_trans,[1 3 2 4]),...
                              [nwaves*numel(spin_system.control.operators) ...
                               numel(spin_system.control.operators)*nwaves]);
            
            % Reorder transformed Hessian
            hess_in_basis(:,:,n)=hess_reorder(hess_trans,nwaves,...
                                              numel(spin_system.control.operators));
                                          
        end
        
        % Write the result out
        hess=hess_in_basis;
        
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,waveform)
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be an array of real numbers.');
end
if size(waveform,1)~=numel(spin_system.control.operators)
    error('the number of rows in waveform must be equal to the number of controls.');
end
if ~isempty(spin_system.control.basis)
    if size(spin_system.control.basis,2)~=spin_system.control.pulse_ntpts
        error(['the number of columns in control.basis must be '...
               int2str(spin_system.control.pulse_ntpts)]);
    end
end
end

% The only "intuitive" interface is the nipple. After 
% that it's all learned.
%
% Bruce Ediger

