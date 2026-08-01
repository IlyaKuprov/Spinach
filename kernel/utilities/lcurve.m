% L-curve analysis function. Syntax:
%
%             lam_opt=lcurve(lam,err,reg,mode)
%
% Parameters:
%
%       lam - row vector of regularisation parameters, must
%             be positive and in ascending order
%
%       err - row vector of least squares errors, must be
%             positive and increasing with lam
%
%       reg - row vector of regularisation functional values,
%             must be positive and decreasing with lam. This
%             is the regularisation functional itself, not the
%             penalty term of the error functional: when the
%             optimiser reports lam*||L*x||^2, divide it by lam
%             once before calling this function
%
%      mode - 'log' for logarithmic coordinates and 'linear'
%             for linear ones; 'log' is recommended
%
% Outputs:
%
%   lam_opt - the regularisation parameter at the point
%             of the maximum curvature of the L-curve
%
% Notes: the corner is the point of the greatest curvature, which for a
%        smooth asymmetric bend is not the intersection of the asymptotes;
%        the criterion locates the regularisation parameter to within a
%        factor of a few and is not convergent in the zero noise limit
%        (Vogel, SIAM J. Numer. Anal. 34, 1996). Treat the answer as an
%        order of magnitude estimate and inspect the plotted curve.
%
%        The curvature maximum must fall inside the sampled interval; if
%        it falls on either end, the corner is outside the range and this
%        function refuses to return the endpoint as an answer.
%
%        This function requires the Curve Fitting Toolbox.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lcurve.m>

function lam_opt=lcurve(lam,err,reg,mode)

% Check consistency
grumble(lam,err,reg,mode);

% Warn about a sweep that is not a trade-off curve, but do not bomb out
if any(diff(err)<=0)||any(diff(reg)>=0)
    warning('err should increase and reg should decrease with lam, inspect the sweep.');
end

% Move to logarithmic coordinates
log_lam=log10(lam); log_err=log10(err); log_reg=log10(reg);

% Resample using quintic spline
sp_err=spapi(optknt(log_lam,5),log_lam,log_err);
sp_reg=spapi(optknt(log_lam,5),log_lam,log_reg);
log_err=fnval(linspace(min(log_lam),max(log_lam),1000),sp_err);
log_reg=fnval(linspace(min(log_lam),max(log_lam),1000),sp_reg);
log_lam=linspace(min(log_lam),max(log_lam),1000); 

% Return to linear coordinates
err=10.^log_err; reg=10.^log_reg; lam=10.^log_lam;

% Plot the L-curve
subplot(1,2,1); plot(err,reg,'b-'); 
hold on; axis tight; kgrid;
kxlabel('least squares error');
kylabel('regularisation error');

% Get the derivatives
switch mode
    
    case 'log'
        
        % Derivatives in logarithmic coordinates
        xp=fdvec(log_err,5,1); xpp=fdvec(log_err,5,2);
        yp=fdvec(log_reg,5,1); ypp=fdvec(log_reg,5,2);
        
        % Plot in logarithmic coordinates
        set(gca,'xscale','log'); set(gca,'yscale','log');
        
    case 'linear'
        
        % Derivatives in linear coordinates
        xp=fdvec(err,5,1); xpp=fdvec(err,5,2);
        yp=fdvec(reg,5,1); ypp=fdvec(reg,5,2);
        
    otherwise
        
        % Complain and bomb out
        error('unknown differentiation mode.');
        
end

% Get the signed curvature
kappa=(xp.*ypp-yp.*xpp)./((xp.^2+yp.^2).^(3/2));

% Plot the curvature
subplot(1,2,2); plot(lam,kappa);
set(gca,'xscale','log'); hold on;
kxlabel('regularisation parameter');
kylabel('L-curve curvature'); 
axis tight; kgrid;

% Find the maximum curvature away from the unreliable stencil ends
margin=ceil(numel(kappa)/40);
[~,index]=max(kappa((1+margin):(end-margin))); index=index+margin;

% A corner at the edge means it is outside the sampled interval
if (index<=(2*margin))||(index>=(numel(kappa)-2*margin))
    error('L-curve corner is outside the regularisation parameter range, widen it.');
end

% Report the optimum point
lam_opt=lam(index);
subplot(1,2,1); plot(err(index),reg(index),'ro');
subplot(1,2,2); plot(lam(index),kappa(index),'ro');

end

% Consistency enforcement
function grumble(lam,err,reg,mode)
if (~isnumeric(lam))||(~isreal(lam))||(any(~isfinite(lam)))||...
   (size(lam,1)~=1)||any(lam<=0)
    error('lam must be a row vector of positive real numbers.');
end
if (~isnumeric(err))||(~isreal(err))||(any(~isfinite(err)))||...
   (size(err,1)~=1)||any(err<=0)
    error('err must be a row vector of positive real numbers.');
end
if (~isnumeric(reg))||(~isreal(reg))||(any(~isfinite(reg)))||...
   (size(reg,1)~=1)||any(reg<=0)
    error('reg must be a row vector of positive real numbers.');
end
if (numel(lam)~=numel(err))||(numel(err)~=numel(reg))
    error('lam, err and reg must have the same number of elements.');
end
if numel(lam)<6
    error('at least six regularisation parameter values are required.');
end
if any(diff(lam)<=0)
    error('lam must be in ascending order.');
end
if (~ischar(mode))||(~ismember(mode,{'log','linear'}))
    error('mode must be either ''log'' or ''linear''.');
end
end

% I don't like ass kissers, flag wavers or team players. I like people who
% buck the system. Individualists. I often warn people: "Somewhere along
% the way, someone is going to tell you, 'There is no "I" in team.' What
% you should tell them is, 'Maybe not. But there is an "I" in independence,
% individuality and integrity.'" Avoid teams at all cost. Keep your circle
% small. Never join a group that has a name. If they say, "We're the So-
% and-Sos," take a walk. And if, somehow, you must join, if it's unavoid-
% able, such as a union or a trade association, go ahead and join. But don't
% participate; it will be your death. And if they tell you you're not a te-
% am player, congratulate them on being observant.
%
% George Carlin

