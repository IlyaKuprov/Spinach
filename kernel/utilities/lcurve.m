% L-curve analysis function. Syntax:
%
%             lambda_opt=lcurve(lambda,err,reg,mode)
%
% Parameters:
%
%       lam - row vector of regularisation parameters
%
%       err - row vector of least squares errors
%
%       reg - row vector of regularisation functional values
%
%      mode - 'log' for logarithmic coordinates and 'linear'
%             for linear ones; 'log' is recommended
%
% The function returns the regularisation parameter at the point 
% of the maximum curvature of the L-curve.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Lcurve.m>

function lam_opt=lcurve(lam,err,reg,mode)

% Check consistency
grumble(lam,err,reg,mode);

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

% Find the optimum point
[~,index]=max(kappa); lam_opt=lam(index);
subplot(1,2,1); plot(err(index),reg(index),'ro');
subplot(1,2,2); plot(lam(index),kappa(index),'ro');

end

% Consistency enforcement
function grumble(lam,err,reg,mode)
if (~isnumeric(lam))||(~isreal(lam))||...
   (any(~isfinite(lam)))||(size(lam,1)~=1)
    error('lam must be a row vector of real numbers.');
end
if (~isnumeric(err))||(~isreal(err))||...
   (any(~isfinite(err)))||(size(err,1)~=1)
    error('err must be a row vector of real numbers.');
end
if (~isnumeric(reg))||(~isreal(reg))||...
   (any(~isfinite(reg)))||(size(reg,1)~=1)
    error('reg must be a row vector of real numbers.');
end
if (numel(lam)~=numel(err))||(numel(err)~=numel(reg))
    error('lam, err and reg must have the same number of elements.');
end
if ~ischar(mode)
    error('mode must be a character string.');
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

