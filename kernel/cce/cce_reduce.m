% Deterministic multiplicative or additive cluster reduction. Syntax:
%
%                 reduced=cce_reduce(raw,clusters,parameters)
%
% raw is a vector of records with cluster_id (0 for central-only), value
% (ntime by nout), abs_error_bound and valid of matching size. Every
% family member appears exactly once; gather order is immaterial. One
% configuration/state is reduced per call, before ensemble averaging.
% Error bounds are absolute and must include propagation uncertainty.
% Invalid input masks propagate as raw_invalid in the reduced result.
%
% parameters.method selects additive noise or multiplicative dynamics;
% reducer is powers (default) or direct. numerics.atol and rtol specify
% nonnegative output budgets (not both zero); zero_policy is mask or
% error (default mask). No epsilon replacement or magnitude clipping.
%
% Output value/valid/reason/abs_error_bound correspond to maximum order;
% by_order and matching diagnostics have an extra final dimension for
% orders 0:max_order. reason codes are ok, raw_invalid, baseline_zero,
% denominator_unresolved, uncertain_zero, arithmetic_range, error_budget.
% Powers use deterministic compensated log-amplitude and phase sums;
% exact zeros with strictly positive powers are valid if their supplied
% uncertainty is zero and every denominator is resolved. Complete bath
% exact-reference bypass belongs to cce.m, not this algebraic reducer.
%
% Inputs:
%   raw - records with cluster_id, time-by-observable value, logical valid,
%   and nonnegative abs_error_bound.
%   clusters - ordered subset-closed records from cce_clusters.
%   parameters - method, numerical tolerances; optional reducer=powers|direct
%   (default powers), numerics.zero_policy=mask|error (default mask).
%
% Outputs:
%   reduced - values, validity/error/reason arrays, by-order results,
%   cluster/cumulant and conditioning diagnostics.
%
% Multiplicative convention: F(empty)=1; F(C)=L(C)/L(empty),
% K(C)=F(C)/product(K(S), nonempty proper S subset C), and
% L_approx=L(empty)*product(K(C), retained nonempty C).
% For noise, K(C)=L(C)-L(empty)-sum(K(S), proper nonempty S),
% and L_approx=L(empty)+sum(K(C)). Masks remain explicit;
% unresolved denominators are never replaced by epsilon.
%
% See also cce_clusters, cce.
% Reference: https://doi.org/10.1103/PhysRevB.78.085315
%
% ilya.kuprov@weizmann.ac.il

function reduced=cce_reduce(raw,clusters,parameters)

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate indexed raw records and the closed cluster family
grumble(raw,clusters,parameters);

% Canonicalise gathered results by stable cluster IDs
[~,baseline_index]=ismember(uint64(0),[raw.cluster_id]);
baseline=raw(baseline_index); nc=numel(clusters);
[~,indices]=ismember([clusters.id],[raw.cluster_id]); records=raw(indices);
shape=size(baseline.value); max_order=max([0 clusters.order]);
by_order=complex(nan([shape max_order+1])); errors=inf([shape max_order+1]);
valid=false([shape max_order+1]); reasons=repmat({'ok'},[shape max_order+1]);
reducer='powers';
if isfield(parameters,'reducer'), reducer=parameters.reducer; end

% Reduce each prefix independently; powers change with the included family
for order=0:max_order

    % Compile exact signed powers of this order's downward-closed family
    active=find([clusters.order]<=order); powers=zeros(1,nc,'int64');
    for k=fliplr(active)
        power=int64(1);
        for j=active(active>k)
            if ismember(k,clusters(j).proper_subsets)
                term=powers(j);
                if (term>0&&power<intmin('int64')+term)||...
                   (term<0&&power>intmax('int64')+term)
                    error('cce_reduce:overflow','signed power overflow.');
                end
                power=power-term;
            end
        end
        powers(k)=power;
    end

    % Evaluate independent times and observable elements in canonical order
    values=complex(nan(shape)); bounds=inf(shape); masks=false(shape);
    why=repmat({'ok'},shape);
    for point=1:numel(values)

        % Gather absolute input uncertainties and validity
        z=complex(zeros(1,nc+1)); eta=zeros(1,nc+1); usable=false(1,nc+1);
        z(1)=baseline.value(point); eta(1)=baseline.abs_error_bound(point);
        usable(1)=baseline.valid(point);
        for k=active
            z(k+1)=records(k).value(point); eta(k+1)=records(k).abs_error_bound(point);
            usable(k+1)=records(k).valid(point);
        end

        % Empty-bath/order-zero result is the unmodified central baseline
        if isempty(active)
            values(point)=z(1); bounds(point)=eta(1); masks(point)=usable(1);
            if ~usable(1)
                why{point}='raw_invalid';
            elseif bounds(point)>parameters.numerics.atol+...
                                 parameters.numerics.rtol*abs(values(point))
                masks(point)=false; why{point}='error_budget';
            end
            continue
        end

        % Additive connected quantities use signed coefficients, not ratios
        if strcmp(parameters.method,'noise')
            terms=[1-sum(double(powers)) double(powers)]; needed=find(terms~=0);
            if any(~usable(needed))
                why{point}='raw_invalid'; continue
            end
            total=complex(0); correction=complex(0); magnitude=0; uncertainty=0;
            for k=needed
                term=terms(k)*z(k); increment=term-correction;
                updated=total+increment; correction=(updated-total)-increment; total=updated;
                magnitude=magnitude+abs(term); uncertainty=uncertainty+abs(terms(k))*eta(k);
            end
            operations=8*numel(needed); unit=eps/2;
            if operations*unit>=1
                why{point}='arithmetic_range'; continue
            end
            values(point)=total;
            bounds(point)=uncertainty+operations*unit/(1-operations*unit)*magnitude;
            if ~isfinite(total)||~isfinite(bounds(point))
                why{point}='arithmetic_range'; continue
            end
            masks(point)=true;
        else

            % Normalisation always requires a resolved nonzero baseline
            if ~usable(1)
                why{point}='raw_invalid'; continue
            elseif abs(z(1))<=eta(1)
                why{point}='baseline_zero'; continue
            end

            % Normalise once and propagate the quotient uncertainty
            g=complex(zeros(1,nc)); ge=zeros(1,nc);
            needed=find(powers~=0);
            if strcmp(reducer,'direct'), needed=active; end
            if any(~usable(needed+1))
                why{point}='raw_invalid'; continue
            end
            unit=eps/2; gamma=6*unit/(1-6*unit);
            for k=needed
                g(k)=z(k+1)/z(1);
                ge(k)=(eta(k+1)+abs(g(k))*eta(1))/(abs(z(1))-eta(1))+gamma*abs(g(k));
            end

            % Use compensated signed log products without phase unwrapping
            if strcmp(reducer,'powers')
                [total,bound,reason]=power_product(z(1),eta(1),g,ge,powers);
                if ~strcmp(reason,'ok')
                    why{point}=reason; continue
                end
            else

                % Evaluate direct irreducible factors in subset order
                factors=complex(ones(1,nc)); factor_errors=zeros(1,nc); reason='ok';
                for k=active
                    subsets=clusters(k).proper_subsets;
                    [denominator,den_error,reason]=power_product(1,0,factors(subsets),...
                                                                factor_errors(subsets),...
                                                                ones(size(subsets),'int64'));
                    if ~strcmp(reason,'ok'), break; end
                    if abs(denominator)<=den_error
                        reason='denominator_unresolved'; break
                    end
                    factors(k)=g(k)/denominator;
                    factor_errors(k)=(ge(k)+abs(factors(k))*den_error)/...
                                     (abs(denominator)-den_error)+gamma*abs(factors(k));
                end
                if strcmp(reason,'ok')
                    [total,bound,reason]=power_product(z(1),eta(1),factors(active),...
                                                      factor_errors(active),...
                                                      ones(size(active),'int64'));
                end
                if ~strcmp(reason,'ok')
                    why{point}=reason; continue
                end

            end
            values(point)=total; bounds(point)=bound; masks(point)=true;

        end

        % Reject estimates exceeding the caller's absolute/relative budget
        if masks(point)&&bounds(point)>parameters.numerics.atol+...
                                       parameters.numerics.rtol*abs(values(point))
            masks(point)=false; why{point}='error_budget';
        end

    end

    % Preserve invalid masks and reason codes, never plausible finite fillers
    values(~masks)=NaN; by_order(:,:,order+1)=values;
    errors(:,:,order+1)=bounds; valid(:,:,order+1)=masks; reasons(:,:,order+1)=why;

end

% Return final order plus every prefix's conditioning diagnostics
reduced=struct('value',by_order(:,:,end),'abs_error_bound',errors(:,:,end),...
               'valid',valid(:,:,end),'reason',{reasons(:,:,end)},...
               'by_order',by_order,'by_order_error',errors,...
               'by_order_valid',valid,'by_order_reason',{reasons},...
               'orders',0:max_order,'reducer',reducer);
if isfield(parameters.numerics,'zero_policy')&&...
   strcmp(parameters.numerics.zero_policy,'error')&&any(~reduced.valid,'all')
    first=find(~reduced.valid,1);
    error('cce_reduce:invalid','unresolved reduction: %s.',reduced.reason{first});
end

end

% Repeated stable signed product with data-dependent rounding/error bounds
function [value,bound,reason]=power_product(base,base_error,z,eta,powers)

% Allocate compensated amplitude, phase and uncertainty accumulators
value=NaN; bound=Inf; reason='ok'; indices=find(powers~=0);
amplitude=log(abs(base)); phase=angle(base); ac=0; pc=0;
condition=base_error/(abs(base)-base_error); exact_zero=false;
scale=abs(amplitude)+abs(phase); count=0;

% Inspect all denominators before accepting any exact numerator zero
for k=indices

    % Refuse nonfinite intermediate quotients before accepting exact zeros
    if ~isfinite(z(k))||~isfinite(eta(k))
        reason='arithmetic_range'; return
    end
    if powers(k)<0&&abs(z(k))<=eta(k)
        reason='denominator_unresolved'; return
    elseif abs(z(k))<=eta(k)
        if z(k)==0&&eta(k)==0&&powers(k)>0
            exact_zero=true; continue
        else
            reason='uncertain_zero'; return
        end
    end
    if abs(powers(k))>flintmax
        reason='arithmetic_range'; return
    end

    % Integer powers preserve complex phase without logarithm branch fixes
    p=double(powers(k)); a=p*log(abs(z(k))); b=p*angle(z(k));
    increment=a-ac; updated=amplitude+increment;
    ac=(updated-amplitude)-increment; amplitude=updated;
    increment=b-pc; updated=phase+increment;
    pc=(updated-phase)-increment; phase=updated;
    condition=condition+abs(p)*eta(k)/(abs(z(k))-eta(k));
    scale=scale+abs(a)+abs(b); count=count+1;
end

% Exact positive zeros need no logarithm or epsilon replacement
if exact_zero
    value=0; bound=0; return
end

% Convert conditioning and rounding bounds before exponentiating
operations=24*count+8; unit=eps/2;
if operations*unit>=1||~isfinite(amplitude)||~isfinite(phase)
    reason='arithmetic_range'; return
end
condition=condition+operations*unit/(1-operations*unit)*(1+scale);
if amplitude>log(realmax)||amplitude<log(realmin*eps)||~isfinite(condition)
    reason='arithmetic_range'; return
end
value=exp(amplitude)*complex(cos(phase),sin(phase));
if condition==0
    bound=0;
else
    log_error=condition+log(-expm1(-condition));
    if amplitude+log_error>log(realmax)
        reason='arithmetic_range'; value=NaN; return
    end
    bound=exp(amplitude+log_error);
end

end

% Validate raw gather indexing, scalar budgets and the subset order
function grumble(raw,clusters,parameters)
if ~isstruct(raw)||isempty(raw)||~isvector(raw)||...
   ~all(isfield(raw,{'cluster_id','value','abs_error_bound','valid'}))||...
   ~isstruct(clusters)||(~isempty(clusters)&&...
   ~all(isfield(clusters,{'id','order','proper_subsets','mobius_power'})))
    error('cce_reduce:inputs','indexed raw records and cluster metadata required.');
end
ids=[raw.cluster_id];
if numel(ids)~=numel(raw)||numel(unique(ids))~=numel(ids)||...
   ~all(isfinite(ids))||any(ids<0)||any(ids~=fix(ids))||sum(ids==0)~=1||...
   numel(raw)~=numel(clusters)+1||~all(ismember([clusters.id],ids))
    error('cce_reduce:ids','raw requires one baseline and exactly one record per cluster.');
end
shape=size(raw(1).value);
for k=1:numel(raw)
    r=raw(k);
    if ~isnumeric(r.value)||~ismatrix(r.value)||isempty(r.value)||...
       ~isequal(size(r.value),shape)||~isnumeric(r.abs_error_bound)||...
       ~isreal(r.abs_error_bound)||~isequal(size(r.abs_error_bound),shape)||...
       any(r.abs_error_bound<0,'all')||any(~isfinite(r.abs_error_bound),'all')||...
       ~islogical(r.valid)||~isequal(size(r.valid),shape)||...
       any(~isfinite(r.value(r.valid)))
        error('cce_reduce:raw','raw values, finite absolute errors and logical masks must match.');
    end
end
if ~issorted([clusters.order])||numel(unique([clusters.id]))~=numel(clusters)
    error('cce_reduce:family','cluster family must have unique IDs in increasing order.');
end
for k=1:numel(clusters)
    s=clusters(k).proper_subsets;
    if ~isscalar(clusters(k).order)||~isfinite(clusters(k).order)||...
       clusters(k).order<1||clusters(k).order~=fix(clusters(k).order)||...
       clusters(k).id==0
        error('cce_reduce:family','nonempty clusters require positive order and ID.');
    end
    if any(s<1)||any(s>=k)||any(s~=fix(s))||numel(unique(s))~=numel(s)||...
       any([clusters(s).order]>=clusters(k).order)||...
       ~isa(clusters(k).mobius_power,'int64')||~isscalar(clusters(k).mobius_power)
        error('cce_reduce:subsets','proper subsets must precede their parent.');
    end
end
if ~isstruct(parameters)||~isscalar(parameters)||...
   ~all(isfield(parameters,{'method','numerics'}))||...
   ~ismember(parameters.method,{'cce','gcce','mecce','megcce','noise'})||...
   ~isstruct(parameters.numerics)||~isscalar(parameters.numerics)||...
   ~all(isfield(parameters.numerics,{'atol','rtol'}))
    error('cce_reduce:parameters','method and numerical budgets required.');
end
for field={'atol','rtol'}
    value=parameters.numerics.(field{1});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<0
        error('cce_reduce:budget','absolute and relative budgets must be finite nonnegative scalars.');
    end
end
if parameters.numerics.atol+parameters.numerics.rtol==0
    error('cce_reduce:budget','absolute and relative budgets cannot both be zero.');
end
if isfield(parameters,'reducer')&&~ismember(parameters.reducer,{'powers','direct'})
    error('cce_reduce:reducer','reducer must be powers or direct.');
end
if isfield(parameters.numerics,'zero_policy')&&...
   ~ismember(parameters.numerics.zero_policy,{'mask','error'})
    error('cce_reduce:policy','zero_policy must be mask or error.');
end
end


