% Explicit-window decay diagnostics without an automatic T1 claim. Syntax:
%
%                    analysis=cce_analyse(result,parameters)
%
% result contains times, value, valid. parameters.analysis.fit defaults
% to false. When true, require observable (column), envelope='magnitude'
% or 'real', model='exponential' or 'stretched_exponential', interval=[a b]
% in seconds, weights (one positive inverse-variance weight per time),
% confidence='none' or 'linearized', and confidence_level in (0,1).
% Model A*exp(-(t/T)^beta) has no offset. Fits use weighted least squares
% in the declared envelope, not its logarithm; confidence intervals use
% local Jacobian Gaussian errors and are not a CCE convergence claim.
% Invalid or nondecaying/oscillatory fitting windows are refused with
% warnings in the returned structure, never silently deleted or clipped.
%
% Inputs:
%   result - times, value and valid.
%   parameters.analysis - fit (default false); when enabled: observable
%   column, envelope, model, interval in seconds, positive weights,
%   confidence and confidence_level.
%
% Outputs:
%   analysis - retained data/masks, status, warnings and fit
%   parameters/diagnostics.
%
% See also cce_diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function analysis=cce_analyse(result,parameters)

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate retained raw data and the explicit fitting contract
grumble(result,parameters);
analysis=struct('times',result.times,'value',result.value,'valid',result.valid,...
                'status','not_requested','warnings',{{}},'fit',[]);
if ~isfield(parameters,'analysis')||~isfield(parameters.analysis,'fit')||...
   ~parameters.analysis.fit, return; end
a=parameters.analysis; t=result.times(:); column=a.observable;
y=result.value(:,column); valid=result.valid(:,column);
if strcmp(a.envelope,'magnitude'), y=abs(y); else, y=real(y); end
window=t>=a.interval(1)&t<=a.interval(2); t=t(window); y=y(window);
w=a.weights(window); w=w(:); np=2+strcmp(a.model,'stretched_exponential');

% Refuse invalid or unidentifiable windows without masking evidence
if numel(t)<=np||any(~valid(window))||any(~isfinite(y))||any(y<=0)
    analysis.status='invalid_window'; analysis.warnings={'Invalid, nonpositive or insufficient fitting samples.'};
    return
end
scale=max(y); increases=diff(y)>sqrt(eps)*scale;
if any(increases)||y(end)>=y(1)
    analysis.status='not_decay'; analysis.warnings={'Oscillatory or nondecaying input: declare a valid envelope/window.'};
    return
end

% Fit a dimensionless time coordinate with strictly positive model parameters
time_scale=t(end); x=t/time_scale;
initial=[log(y(1));log(max(eps,-(x(end)-x(1))/log(y(end)/y(1))))];
if np==3, initial(3)=0; end
objective=@(q) sum(w.*(decay(q,x)-y).^2)/sum(w);
options=optimset('Display','off','TolX',1e-10,'TolFun',eps*scale^2,...
                 'MaxFunEvals',5000,'MaxIter',3000);
[q,~,flag]=fminsearch(objective,initial,options); fitted=decay(q,x);
residual=y-fitted; beta=1;
if np==3, beta=exp(q(3)); end
fit=struct('amplitude',exp(q(1)),'time_constant',exp(q(2))*time_scale,...
           'beta',beta,'model',a.model,'interval',a.interval,...
           'envelope',a.envelope,'observable',column,'times',t,...
           'fitted',fitted,'residual',residual,'weights',w,...
           'weighted_rms',sqrt(sum(w.*residual.^2)/sum(w)),...
           'confidence',a.confidence,'confidence_level',a.confidence_level,...
           'confidence_interval',[],'optimizer_exitflag',flag);

% Compute optional local Gaussian confidence intervals in positive coordinates
if strcmp(a.confidence,'linearized')
    power=(x/exp(q(2))).^beta; jac=[fitted fitted.*beta.*power];
    if np==3
        logarithm=zeros(size(x)); nonzero=x>0; logarithm(nonzero)=log(x(nonzero)/exp(q(2)));
        jac(:,3)=-fitted.*power.*beta.*logarithm;
    end
    information=jac'*(w.*jac);
    if rcond(information)<=eps
        analysis.warnings{end+1}='Rank-deficient confidence Jacobian.';
    else
        variance=sum(w.*residual.^2)/(numel(t)-np);
        covariance=variance*(information\eye(np));
        z=sqrt(2)*erfinv(a.confidence_level); radius=z*sqrt(max(0,diag(covariance)));
        intervals=exp([q-radius q+radius]); intervals(2,:)=intervals(2,:)*time_scale;
        fit.confidence_interval=intervals;
    end
end

% Preserve the selected observable and fit assumptions with an explicit status
analysis.fit=fit; analysis.status='fitted';
if flag<=0
    analysis.status='fit_limit'; analysis.warnings{end+1}='Optimizer iteration or evaluation limit reached.';
end

end

% Evaluate the positive-parameter envelope in scaled seconds
function values=decay(q,t)
beta=1;
if numel(q)==3, beta=exp(q(3)); end
values=exp(q(1))*exp(-(t/exp(q(2))).^beta);
end

% Check observations, confidence procedure and fitting intervals
function grumble(result,parameters)
if ~isstruct(result)||~all(isfield(result,{'times','value','valid'}))||...
   ~isnumeric(result.times)||~isreal(result.times)||~isvector(result.times)||...
   any(~isfinite(result.times))||any(result.times<0)||any(diff(result.times)<=0)||...
   ~isnumeric(result.value)||size(result.value,1)~=numel(result.times)||...
   ~islogical(result.valid)||~isequal(size(result.valid),size(result.value))||...
   ~isstruct(parameters)
    error('Result requires increasing nonnegative times, values and matching logical validity.');
end
if ~isfield(parameters,'analysis'), return; end
a=parameters.analysis;
if ~isstruct(a)||~isfield(a,'fit')||~islogical(a.fit)||~isscalar(a.fit)
    error('analysis.fit must be a scalar logical.');
end
if ~a.fit, return; end
if ~all(isfield(a,{'observable','envelope','model','interval','weights',...
                  'confidence','confidence_level'}))||...
   ~isscalar(a.observable)||a.observable~=fix(a.observable)||...
   a.observable<1||a.observable>size(result.value,2)||...
   ~ismember(a.envelope,{'magnitude','real'})||...
   ~ismember(a.model,{'exponential','stretched_exponential'})||...
   ~isnumeric(a.interval)||~isreal(a.interval)||numel(a.interval)~=2||any(~isfinite(a.interval))||...
   a.interval(1)<0||a.interval(2)<=a.interval(1)||...
   ~isnumeric(a.weights)||~isreal(a.weights)||~isvector(a.weights)||...
   numel(a.weights)~=numel(result.times)||any(~isfinite(a.weights))||any(a.weights<=0)||...
   ~ismember(a.confidence,{'none','linearized'})||...
   ~isscalar(a.confidence_level)||~isreal(a.confidence_level)||...
   ~isfinite(a.confidence_level)||a.confidence_level<=0||a.confidence_level>=1
    error('Fitting requires explicit valid observable, envelope, model, window, weights and confidence.');
end
end


