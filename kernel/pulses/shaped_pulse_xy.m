% Shaped pulse function using Cartesian coordinates. Applies a user-
% specified pulse shape on user-specified operators while the rest of
% the drift Liouvillian continues to affect the spin system. Syntax:
%
%  [rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,...
%                               amplitudes,slice_durs,rho,method)
% 
% Parameters:
%
%       drift - the drift Liouvillian, the part of the Liouvillian that
%               should continue running in the background. This should 
%               include the transmitter offset term, if any. 
%
%    controls - a cell array of control operators corresponding to each
%               channel, this may include operators for spatial degrees
%               of freedom, such as gradients and diffusion.
%
%  amplitudes - a cell array of control amplitude vectors in rad/s, one
%               vector per control channel; the elements of each vector
%               correspond to different time points.
%
%  slice_durs - a vector containing the duration of each pulse slice,
%               seconds. For piecewise-constant methods, the number of
%               durations should be equal to the nuber of amplitudes.
%               For piecewise-linear methods, there should be one ele-
%               ment more in the amplitude array.
%             
%         rho - initial state vector or a bookshelf matrix thereof
%
%      method - propagation method and product quadrature: 
%
%                 Krylov algorithm (usually faster for calls with one 
%                 and two outputs):
%
%                    'expv-pwc' - piecewise-constant
%                    'expv-pwl' - 2nd order Lie quadrature
%
%                 Explicit matrix exponentiation (usually faster for 
%                 calls with three outputs):
%
%                    'expm-pwc' - piecewise-constant
%                    'expm-pwl' - 2nd order Lie quadrature
%
%                 Spinach evolution function call (do not choose un-
%                 less you have a specific good reason):
%                 
%                    'evol-pwc' - piecewise-constant
%                    'evol-pwl' - 2nd order Lie quadrature
%
% Outputs:
%
%         rho - state vector for the final state, or a stack thereof
%
%        traj - system trajectory as a [1 x (nsteps+1)] cell array, 
%               the first point is the initial condition
%
%           P - effective pulse propagator (expensive, best avoided)
%
% i.kuprov@soton.ac.uk
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=shaped_pulse_xy.m>

function [rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,...
                                      amplitudes,slice_durs,rho,method)

% Check consistency
grumble(drift,controls,amplitudes,slice_durs,rho,method);

% Start feedback timer
feedback=tic();

% Get trajectory started
if nargout>1

    % Preallocate as a cell array
    traj=cell(1,numel(slice_durs)+1);

    % First element is
    % initial condition
    traj{1}=rho;

end

% Get the propagator started
if nargout>2, P=speye(size(drift)); end

% Move appropriate things to GPU
if ismember('gpu',spin_system.sys.enable)
    drift=gpuArray(drift);
    for n=1:numel(controls)
        controls{n}=gpuArray(controls{n});
    end
    rho=gpuArray(rho);
    if nargout>2, P=gpuArray(P); end
end

% Decide propagation method
switch method(1:4)
            
    % Krylov method
    case 'expv'

        % Loop over pulse slices
        for n=1:numel(slice_durs)

            % Decide the quadrature
            switch method(6:8)

                % Piecewise-constant
                case 'pwc'

                    % Grab the drift
                    slice_oper=drift;

                    % Add controls
                    for k=1:numel(controls)
                        slice_oper=slice_oper+amplitudes{k}(n)*controls{k};
                    end

                    % Apply the evolution slice
                    rho=step(spin_system,slice_oper,rho,slice_durs(n));

                    % Update the propagator
                    if nargout>2
                        P=step(spin_system,slice_oper,P,slice_durs(n));
                        P=clean_up(spin_system,P,spin_system.tols.prop_chop);
                    end

                % 2nd order Lie
                case 'pwl'

                    % Grab left and right edge drifts
                    slice_oper_l=drift; slice_oper_r=drift;

                    % Add controls
                    for k=1:numel(controls)
                        slice_oper_l=slice_oper_l+amplitudes{k}(n)*controls{k};
                        slice_oper_r=slice_oper_r+amplitudes{k}(n+1)*controls{k};
                    end

                    % Apply the evolution slice
                    rho=step(spin_system,{slice_oper_l,slice_oper_r},rho,slice_durs(n));

                    % Update pulse propagator
                    if nargout>2
                        P=step(spin_system,{slice_oper_l,slice_oper_r},P,slice_durs(n));
                        P=clean_up(spin_system,P,spin_system.tols.prop_chop);
                    end

                otherwise

                    % Complain and bomb out
                    error('unknown product quadrature type');

            end

            % Store trajectory element
            if nargout>1, traj{n+1}=rho; end

            % Inform the user
            if (n==numel(slice_durs))||(toc(feedback)>1)

                % Report the progress
                report(spin_system,['shaped pulse slice ' num2str(n) ...
                                    '/' num2str(numel(slice_durs)) '...']);

                % Restart timer
                feedback=tic();

            end

        end

    % Matrix exponentiation
    case 'expm'
                
        % Loop over pulse slices
        for n=1:numel(slice_durs)

            % Decide the quadrature
            switch method(6:8)

                % Piecewise-constant
                case 'pwc'

                    % Grab the drift
                    slice_oper=drift;

                    % Add controls
                    for k=1:numel(controls)
                        slice_oper=slice_oper+amplitudes{k}(n)*controls{k};
                    end
                    
                case 'pwl'

                    % Grab left and right point drifts
                    slice_oper_l=drift; slice_oper_r=drift;

                    % Add controls
                    for k=1:numel(controls)
                        slice_oper_l=slice_oper_l+amplitudes{k}(n)*controls{k};
                        slice_oper_r=slice_oper_r+amplitudes{k}(n+1)*controls{k};
                    end

                    % Build the evolution generator
                    slice_oper=isergen(slice_oper_l,[],slice_oper_r,slice_durs(n));

                otherwise

                    % Complain and bomb out
                    error('unknown product quadrature type');

            end

            % Compute slice propagator
            PS=propagator(spin_system,slice_oper,slice_durs(n));

            % Update state
            rho=PS*rho;

            % Store trajectory element
            if nargout>1, traj{n+1}=rho; end

            % Update pulse propagator
            if nargout>2
                P=clean_up(spin_system,PS*P,spin_system.tols.prop_chop);
            end

            % Inform the user
            if (n==numel(slice_durs))||(toc(feedback)>1)

                % Report the progress
                report(spin_system,['shaped pulse slice ' num2str(n) ...
                                    '/' num2str(numel(slice_durs)) '...']);

                % Restart timer
                feedback=tic();

            end

        end

    % Spinach evolution() call
    case 'evol'
                
        % Loop over pulse slices
        for n=1:numel(slice_durs)

            % Decide the quadrature
            switch method(6:8)

                % Piecewise-constant
                case 'pwc'

                    % Grab the drift
                    slice_oper=drift;

                    % Add controls
                    for k=1:numel(controls)
                        slice_oper=slice_oper+amplitudes{k}(n)*controls{k};
                    end
                    
                case 'pwl'

                    % Grab left and right point drifts
                    slice_oper_l=drift; slice_oper_r=drift;

                    % Add controls
                    for k=1:numel(controls)
                        slice_oper_l=slice_oper_l+amplitudes{k}(n)*controls{k};
                        slice_oper_r=slice_oper_r+amplitudes{k}(n+1)*controls{k};
                    end

                    % Build the evolution generator
                    slice_oper=isergen(slice_oper_l,[],slice_oper_r,slice_durs(n));

                otherwise

                    % Complain and bomb out
                    error('unknown product quadrature type');

            end

            % Evolve the state
            rho=evolution(spin_system,slice_oper,[],rho,slice_durs(n),1,'final');

            % Store trajectory element
            if nargout>1, traj{n+1}=rho; end

            % Update pulse propagator
            if nargout>2
                P=evolution(spin_system,slice_oper,[],P,slice_durs(n),1,'final');
                P=clean_up(spin_system,P,spin_system.tols.prop_chop);
            end

            % Inform the user
            if (n==numel(slice_durs))||(toc(feedback)>1)

                % Report the progress
                report(spin_system,['shaped pulse slice ' num2str(n) ...
                                    '/' num2str(numel(slice_durs)) '...']);

                % Restart timer
                feedback=tic();

            end

        end

    otherwise

        % Complain and bomb out
        error('unknown propagation method.');

end

% Retrieve the outputs from GPU
if ismember('gpu',spin_system.sys.enable)
    rho=gather(rho);
    if nargout>1
        for n=1:numel(traj)
            traj{n}=gather(traj{n});
        end
    end
    if nargout>2, P=gather(P); end
end

end

% Consistency enforcement
function grumble(drift,controls,amplitudes,slice_durs,rho,method)
if (~isnumeric(drift))||(size(drift,1)~=size(drift,2))
    error('drift operator must be a square matrix.');
end
if (~iscell(controls))||isempty(controls)
    error('controls must be a cell array of matrices.');
end
for n=1:numel(controls)
    if ~isnumeric(controls{n})
        error('control operators must be matrices.');
    end
    if (size(controls{n},1)~=size(controls{n},2))||...
       (size(controls{n},1)~=size(drift,1))
        error('control operators must have the same size as drift.');
    end
end
if (~isnumeric(slice_durs))||(~isreal(slice_durs))||...
     any(~isfinite(slice_durs),'all')||any(slice_durs<0,'all')
    error('slice_durs must be a vector of non-negative real numbers.');
end
if (~iscell(amplitudes))||isempty(amplitudes)
    error('amplitudes must be a cell array of vectors.');
end
for n=1:numel(amplitudes)
    if ~isnumeric(amplitudes{n})
        error('all elements of amplitudes array must be vectors.');
    end
    if strcmp(method(6:8),'pwl')&&(numel(amplitudes{n})~=numel(slice_durs)+1)
        error(['slice_durs must have one fewer element than amplitude '   ...
               'tables for 2nd order Lie product quadrature.']);
    end
    if strcmp(method(6:8),'pwc')&&(numel(amplitudes{n})~=numel(slice_durs))
        error(['slice_durs must have the same number of elements as amplitude '   ...
               'tables for piecewise-constant product quadrature.']);
    end
end
if size(drift,2)~=size(rho,1)
    error('the dimensions of drift and rho must agree.');
end
end

% When IK proposed the fibre etching technique featured in his 2004
% JMR paper (http://dx.doi.org/10.1016/j.jmr.2004.08.017), he could
% see terror in his supervisors's eyes -- Peter Hore reasonably tho-
% ught that the notoriously eccentric Russian student could not pos-
% sibly be trusted with boiling hydrofluoric acid in a high-power 
% laser lab. The Chemistry Department safety officer held a similar
% view. It is not entirely clear how IK got hold of several millili-
% tres of concentrated HF and a heat gun on a Saturday night in the
% PTCL Teaching Lab, but the photographs of the resulting fibre tip
% were left on Peter's table on Monday. The paper was accepted by 
% JMR without revisions.

