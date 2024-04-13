% Penalty terms for the Optimal Control module. Returns the penalty
% function and its gradient for the waveform, which should be sup-
% plied as row vector or a horizontal stack thereof. Syntax:
%
%        [pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)
% 
% Parameters:
%
%    wf                          -  control sequence waveform
%
%    type='none'                 -  no waveform penalty.
%
%    type='NS'                   -  norm square, designed to favour
%                                   low-power waveforms over high-
%                                   power ones.
%
%    type='DNS'                  -  derivative norm square, desig-
%                                   ned to favour smooth waveforms
%                                   over jagged ones.
%
%    type='SNS'                  -  spillout norm square, NS appli-
%                                   ed to the part of the waveform
%                                   with values outside the floor
%                                   and ceiling bounds.
%
%    type='SNSA'                 -  SNS applied after a transform to
%                                   amplitude-phase representation
%                                   Penalises amplitude values outs-
%                                   ide the ceiling bound. Requires
%                                   even number of control channels
%                                   with waveform rows ordered as:
%                                   [Xa Ya Xb Yb ... Xn Yn]
%
%    fb                          -  floor bound, the lower bound on
%                                   waveform used in the SNS penal-
%                                   ty function.
%                                   
%    cb                          -  ceiling bound, the upper bound
%                                   on waveform used in the SNS and
%                                   SNSA penalty functions.
%
% Outputs:
%
%    pen_term                    -  value of the penalty term
%
%    pen_grad                    -  gradient of the penalty term with
%                                   respect to the waveform vector
%
%    pen_hess                    -  Hessian of the penalty term with
%                                   respect to the waveform vector 
%
% The waveforms on different channels are assumed to be stored in the 
% rows of the input array. The Hessian elements correspond to the ele-
% ments of the waveform array ordered as:
%
%                  [X1 Y1 Z1 X2 Y2 Z2 ... Xn Yn Zn]
%
% where X,Y,Z are different control channels and the index enumerates
% the time discretization points. Gradient dimensions and element or-
% der are the same as the input waveform dimensions and element order.
% 
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=penalty.m>

function [pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)

% Check consistency
grumble(wf,type,fb,cb);

% Preallocate the results
if nargout>0, pen_term=0; end
if nargout>1, pen_grad=zeros(size(wf)); end
if nargout>2, pen_hess=zeros(numel(wf),numel(wf)); end

% Decide penalty type
switch type
    
    case 'none'
        
        % Nothing to do
    
    case 'NS'    % Norm square
        
        % Compute the penalty
        if nargout>0
            pen_term=sum(sum(wf.^2));
            pen_term=pen_term/size(wf,2);
        end
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*wf;
            pen_grad=pen_grad/size(wf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*speye(numel(wf));
            pen_hess=pen_hess/size(wf,2);
        end
    
    case 'DNS'   % Derivative norm square
        
        % Five-point second derivative matrix
        D=fdmat(size(wf,2),5,2,'wall');
        
        % Compute the penalty
        if nargout>0
            dwf=wf*D';
            pen_term=sum(sum(dwf.^2));
            pen_term=pen_term/size(dwf,2);
        end
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*dwf*D;
            pen_grad=pen_grad/size(dwf,2);
        end
            
        % Compute the Hessian
        if nargout>2
            pen_hess=2*kron(D'*D,speye(size(wf,1)));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'SNS'   % Spillout norm square
        
        % Build ceiling hole inventory
        ch_map=(wf>cb); ch_actual=wf.*ch_map; ch_wanted=cb.*ch_map;
        
        % Build floor hole inventory
        fh_map=(wf<fb); fh_actual=wf.*fh_map; fh_wanted=fb.*fh_map;
        
        % Compute the penalty
        pen_term=pen_term+sum(sum((ch_actual-ch_wanted).^2))+...
                          sum(sum((fh_actual-fh_wanted).^2));
        pen_term=pen_term/size(wf,2);
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*(ch_actual-ch_wanted)+...
                     2*(fh_actual-fh_wanted);
            pen_grad=pen_grad/size(wf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*ch_map+2*fh_map;
            pen_hess=diag(pen_hess(:));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'SNSA'  % Spillout norm square on amplitudes
        
        % Preallocate amplitude and phase vectors
        amp=zeros(size(wf,1)/2,size(wf,2));
        phi=zeros(size(wf,1)/2,size(wf,2));
        
        % Calculate amplitude and phase
        for n=1:size(wf,1)/2
            [amp(n,:),phi(n,:)]=cartesian2polar(wf(2*n-1,:),wf(2*n,:));
        end
        
        % Build the amplitude hole inventory
        ch_map=(amp>cb); ch_wanted=cb.*ch_map;
        
        % Preallocated new Cartesian bounds
        b=zeros(size(wf));
        
        % Transform to Cartesian representation
        for n=1:size(wf,1)/2
            [b(2*n-1,:),b(2*n,:)]=polar2cartesian(ch_wanted(n,:),phi(n,:));
        end
        
        % New Cartesian bounds from amplitude bounds
        cb=b; cb(cb<=0)=max(max(wf))+1;
        fb=b; fb(fb>=0)=min(min(wf))-1;
        
        % Call the SNS penalty with new bounds
        if nargout==1
            [pen_term]=penalty(wf,'SNS',fb,cb);
        elseif nargout==2
            [pen_term,pen_grad]=penalty(wf,'SNS',fb,cb);
        elseif nargout==3
            [pen_term,pen_grad,pen_hess]=penalty(wf,'SNS',fb,cb);
        end
       
    otherwise
        
        % Complain and bomb out
        error('unknown penalty function type.');
        
end

end

% Consistency enforcement
function grumble(wf,type,fb,cb)
if ~isnumeric(wf)||(~isreal(wf))
    error('wf must be a real numeric array.');
end
if ~isnumeric(fb)||(~isreal(fb))
    error('fb must be a real numeric array.');
end
if ~isnumeric(cb)||(~isreal(cb))
    error('cb must be a real numeric array.');
end
if ~(isscalar(fb)||(isnumeric(fb)&&isequal(size(fb),size(wf))))
    error('floor bound must be a scalar or a numeric array.');
end
if ~(isscalar(cb)||(isnumeric(cb)&&isequal(size(cb),size(wf))))
    error('ceiling bound must be a scalar or a numeric array.');
end
if any(cb<fb)
    error('ceiling bound should be above the floor bound.');
end
if ~ischar(type)
    error('type must be a character string.');
end
if ismember(type,{'SNSA'})&&logical(mod(size(wf,1),2))
    error('waveform must have an even number of rows.')
end
end

% I would never die for my beliefs because I might be wrong.
%
% Bertrand Russell

