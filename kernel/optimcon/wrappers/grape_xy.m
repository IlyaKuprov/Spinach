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
% ilya.kuprov@weizmann.ac.il
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
        grad=tensorprod(spin_system.control.basis,grad,2,2);
        grad=permute(grad,[2 1 3]);
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
        
        % Transform gradients
        grad=tensorprod(spin_system.control.basis,grad,2,2);
        grad=permute(grad,[2 1 3]);
        
        % Preallocate transformed Hessians
        hess_in_basis=zeros([nwaves*numel(spin_system.control.operators) ...
                             nwaves*numel(spin_system.control.operators) ...
                             npenterms+1]);

        % Transform Hessians
        for n=1:size(hess,3)

            % Reorder Hessian
            hess_re=hess_reorder(hess(:,:,n),numel(spin_system.control.operators),...
                                             spin_system.control.pulse_nsteps);

            % Fold up the block matrix into a 4D tensor
            hess_re=reshape(hess_re,[spin_system.control.pulse_nsteps ...
                                     numel(spin_system.control.operators) ...
                                     spin_system.control.pulse_nsteps ...
                                     numel(spin_system.control.operators)]);

            % Transform the Hessian
            hess_re=tensorprod(spin_system.control.basis,hess_re,2,1);
            hess_re=tensorprod(spin_system.control.basis,hess_re,2,3);
            hess_re=permute(hess_re,[2 3 1 4]);

            % Unfold the 4D tensor into block matrix
            hess_re=reshape(hess_re,[nwaves*numel(spin_system.control.operators) ...
                                     nwaves*numel(spin_system.control.operators)]);

            % Reorder transformed Hessian
            hess_in_basis(:,:,n)=hess_reorder(hess_re,nwaves,...
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

