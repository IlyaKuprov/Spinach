% Finds a local maximum of a function of several variables using Newton
% and quasi-Newton algorithms. Syntax:
%
%          [x,data]=fmaxnewton(spin_system,cost_function,guess)
%
% Parameters:
%
%    spin_system   - Spinach data object that has been
%                    through optimcon.m function
%
%    cost_function - a function handle that takes the input
%                    the size of guess
%
%    guess         - the initial point of the optimisation
%
% Outputs:
%
%    x               - the final point of the optimisation
%
%    data.count.iter - iteration counter
%
%    data.count.fx   - function evaluation counter
%
%    data.count.gfx  - gradient evaluation counter
%
%    data.count.hfx  - Hessian evaluation counter
%
%    data.count.rfo  - RFO iteration counter
%
%    data.x_shape    - output of size(guess)
%
%    data.*          - further fields may be set by the 
%                      objective functon
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fmaxnewton.m>

function [x,data]=fmaxnewton(spin_system,cost_function,guess)

% Check consistency
grumble(spin_system,cost_function,guess);

% Initialise counters
data.count.iter=0; data.count.fx=0;
data.count.gfx=0;  data.count.hfx=0;
data.count.rfo=0;        

% Look for checkpoint file
if isfield(spin_system.control,'checkpoint')&&...
   exist([spin_system.sys.scratch filesep ...
          spin_system.control.checkpoint],'file')

    % Read the initial guess from checkpoint file in global scratch
    report(spin_system,'WARNING: optimisation restarted from the checkpoint file');
    load([spin_system.sys.scratch filesep ...
          spin_system.control.checkpoint],'x'); 
    if numel(x)~=numel(guess)
        error('waveform size mismatch between guess and checkpoint.');
    end
    data.x_shape=size(guess);
    
else
    
    % Stretch the guess supplied by user
    data.x_shape=size(guess); x=guess(:);
    
end

% Prepare the freeze mask
if ~isempty(spin_system.control.freeze)

    % Validate freeze mask
    if ~all(size(spin_system.control.freeze)==data.x_shape)
        error(['freeze mask must have dimensions ' int2str(data.x_shape)]);
    end

    % Stretch freeze mask to match x
    frozen=spin_system.control.freeze(:);

else

    % Do not freeze anything
    frozen=false(size(x));

end

% Print the header
header(spin_system);

% Get the video writer going (use VLC to play)
if isfield(spin_system.control,'video_file')

    % This format is available on all platforms
    VW=VideoWriter(spin_system.control.video_file,...
                   'Motion JPEG AVI'); open(VW);

end

% If zero iterations, still display graphics
if spin_system.control.max_iter==0
    [~,~]=objeval(x,cost_function,data,spin_system);
end

% Start the iteration loop
for n=1:spin_system.control.max_iter

    % Default exit flag
    exitflag=0;
    
    % Update iteration counter
    data.count.iter=data.count.iter+1;
    
    % Get the search direction
    switch spin_system.control.method
        
        case 'lbfgs'
            
            if n==1
                
                % Get objective and gradient
                [data,fx,g]=objeval(x,cost_function,data,spin_system);

                % Catch unreasonably small initial fidelities and gradients
                if (abs(data.fx_sep_pen(1))<1e-6)||(norm(g(~frozen),2)<1e-6)
                    error('fidelity or gradient too small at iter 1, find a better guess.');
                end
                
                % Start history arrays
                old_x=x(~frozen); dx_hist=[]; 
                old_g=g(~frozen); dg_hist=[]; 
                
                % Take a conservative step at first iteration
                dir=g.*(~frozen); dir=0.01*dir/max(abs(dir));
                
            else
                
                % Update the history of dx and dg
                dx_hist=[x(~frozen)-old_x dx_hist]; old_x=x(~frozen); %#ok<AGROW>
                dg_hist=[g(~frozen)-old_g dg_hist]; old_g=g(~frozen); %#ok<AGROW>
                
                % Truncate the history
                if size(dx_hist,2)>spin_system.control.n_grads
                    dx_hist=dx_hist(:,1:spin_system.control.n_grads);
                end
                if size(dg_hist,2)>spin_system.control.n_grads
                    dg_hist=dg_hist(:,1:spin_system.control.n_grads);
                end
                
                % Get ascent direction
                dir=zeros(size(frozen));
                dir(~frozen)=lbfgs(dx_hist,dg_hist,g(~frozen));
                
            end
            
        case {'newton','goodwin'}
            
            % Get objective, gradient, and Hessian
            [data,fx,g,H]=objeval(x,cost_function,data,spin_system);
            
            % Tidy up the inputs
            H=real(H+H')/2; g=real(g);
            
            % Apply freeze mask
            H=H(~frozen,~frozen);
            
            % Regularise the Hessian
            [H,data]=hessreg(spin_system,-H,g(~frozen),data);
            
            % Get the search direction
            dir=zeros(size(frozen));
            dir(~frozen)=H\g(~frozen);
            
    end

    % Store the reference point
    g_ref=g(~frozen); dir_ref=dir(~frozen);

    % If line search would be worthwhile, get a bracket [A B] of acceptable points
    [A,B,alpha,fx_new,g_new,next_act,data]=bracketing(cost_function,1,dir,x,fx,...
                                                      g.*(~frozen),data,spin_system);

    % Run sectioning if necessary
    if strcmp(next_act,'sectioning')
    
        % Find an acceptable point within the [A B] bracket    
       [alpha,fx,g,exitflag,data]=sectioning(cost_function,A,B,x,fx,g.*(~frozen),...
                                             dir,data,spin_system);
                                 
    else

        % History update
        fx=fx_new; g=g_new;

    end

    % Report reference point diagnostics to user
    itrep(spin_system,fx,g_ref,dir_ref,alpha,data);
                                         
    % If all good
    if exitflag~=-2
        
        % Update x
        x=x+alpha*dir; 

        % Save checkpoint
        if isfield(spin_system.control,'checkpoint')
            save([spin_system.sys.scratch filesep ... 
                  spin_system.control.checkpoint],'x','-v7.3','-nocompression');
        end

    end
    
    % If all good
    if exitflag==0

        % Check termination conditions
        if norm(alpha*dir,1)<spin_system.control.tol_x

            % Step size
            exitflag=2;

        elseif norm(g(~frozen),2)<spin_system.control.tol_g

            % Grad norm
            exitflag=1;

        end

    end

    % Grab the frame
    if isfield(spin_system.control,'video_file')
        writeVideo(VW,getframe(gcf));
    end
       
    % Exit if necessary
    if exitflag, break; end
   
end

% When no iterations were taken
if ~exist('exitflag','var'), exitflag=0; end

% Fold back the waveform
x=reshape(x,data.x_shape);

% Print the footer
footer(spin_system,exitflag,data);

% Shut down the video writer
if isfield(spin_system.control,'video_file'), close(VW); end

end

% Header printing function
function header(spin_system)
report(spin_system,'==============================================================================================');
report(spin_system,'Iter  #f   #g   #H   #R    fidelity      penalties     total        alpha     |grad|     DGA  ');
report(spin_system,'----------------------------------------------------------------------------------------------');
end

% Footer printing function
function footer(spin_system,exitflag,data)
report(spin_system,'----------------------------------------------------------------------------------------------');
switch(spin_system.control.method)
    case 'lbfgs',  data.algorithm='LBFGS method';
    case 'newton', data.algorithm='Newton-Raphson method';
    case 'goodwin', data.algorithm='Newton-Raphson method with Goodwin acceleration';
end
switch exitflag
    case  1, message='norm(gradient,2) < tol_gfx';
    case  2, message='norm(step,1) < tol_x';
    case  0, message='number of iterations exceeded';
    case -2, message='line search found no maximum';
end
report(spin_system,['    Algorithm Used     : ' data.algorithm]);
report(spin_system,['    Exit message       : ' message]);
report(spin_system,['    Iterations         : ' int2str(data.count.iter)]);
report(spin_system,['    Function Count     : ' int2str(data.count.fx)]);
report(spin_system,['    Gradient Count     : ' int2str(data.count.gfx)]);
report(spin_system,['    Hessian Count      : ' int2str(data.count.hfx)]);
report(spin_system,'========================================================================================');
end

% Iteration report function
function itrep(spin_system,fx,g,dir,alpha,data)

% Performance figures
fid=data.fx_sep_pen(1);
pens=sum(data.fx_sep_pen(2:end));

% Angle between search direction and gradient
dga=abs(acosd((g'*dir)/(norm(g,2)*norm(dir,2))));

% Print iteration data
report(spin_system,[pad(num2str(data.count.iter,'%4.0f'),6),...
                    pad(num2str(data.count.fx,'%4.0f'),5),...
                    pad(num2str(data.count.gfx,'%4.0f'),5),...
                    pad(num2str(data.count.hfx,'%4.0f'),5),...
                    pad(num2str(data.count.rfo,'%4.0f'),5),...
                    pad(num2str(fid,'%+9.6f'),11),'   '...
                    pad(num2str(pens,'%+9.6f'),11),'   '...
                    pad(num2str(fx,'%+9.6f'),11),'   '...
                    pad(num2str(alpha,'%4.0e'),9),...
                    pad(num2str(norm(g(:),2),'%0.2e'),12),...
                    pad(num2str(dga,'%9.1f'),10)]);
                
end

% Consistency enforcement
function grumble(spin_system,cost_function,guess) 
if ~isa(cost_function,'function_handle')
    error('cost_function must be a function handle.');
end
if ~isfield(spin_system,'control')
    error('control field is missing from spin_system structure, run optimcon() first.');
end
if (~isnumeric(guess))||(~isreal(guess))
    error('guess must be an array of real numbers.');
end
switch spin_system.control.integrator
    case 'rectangle'
        if isempty(spin_system.control.basis)
            if size(guess,2)~=spin_system.control.pulse_nsteps
                error('the number of columns in guess must be equal to the number of time steps.');
            end
        else
            if size(spin_system.control.basis,2)~=spin_system.control.pulse_nsteps
                error('the number of columns in waveform basis must be equal to the number of time steps.');
            end
            if size(guess,2)~=size(spin_system.control.basis,1)
                error('the number of columns in guess must be equal to the number of basis functions.');
            end
        end 
    case 'trapezium'
        if isempty(spin_system.control.basis)
            if size(guess,2)~=(spin_system.control.pulse_nsteps+1)
                error('the number of columns in guess must be (number of time steps)+1.');
            end
        else
            if size(spin_system.control.basis,2)~=(spin_system.control.pulse_nsteps+1)
                error('the number of columns in waveform basis must be (number of time steps)+1.');
            end
            if size(guess,2)~=size(spin_system.control.basis,1)
                error('the number of columns in guess must be equal to the number of basis functions.');
            end
        end
    otherwise
        error('unknown time propagation algorithm.');
end
end

% There is a sacred horror about everything grand. It is easy to 
% admire mediocrity and hills; but whatever is too lofty, a geni-
% us as well as a mountain, an assembly as well as a masterpiece,
% seen too near, is appalling... People have a strange feeling of
% aversion to anything grand. They see abysses, they do not see
% sublimity; they see the monster, they do not see the prodigy.
%
% Victor Hugo - Ninety-three
