% Deterministic multiplicative or additive cluster reduction. Syntax:
%
%                 reduced=cce_reduce(raw,clusters,parameters)
%
% raw is a vector of records with cluster_id (0 for observed-spin-only), value
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
% exact-reference bypass belongs to the native ensemble reduction. Native
% cce.bath_spins is the graph-selected domain; complete_record indexes the
% corresponding raw record, or is empty if the family is truncated.
%
% Parameters:
%
%   raw - records with cluster_id, time-by-observable value, logical valid,
%   and nonnegative abs_error_bound.
%   clusters - ordered subset-closed records from cce_clusters.
%   parameters - method, numerical tolerances; optional reducer=powers|direct
%   (default powers), numerics.zero_policy=mask|error (default mask).
%
% Outputs:
%
%   reduced - values, validity/error/reason arrays, by-order results,
%   cluster/cumulant and conditioning diagnostics.
%
% Multiplicative convention: F(empty)=1; F(C)=E(C)/E(empty),
% K(C)=F(C)/product(K(S), nonempty proper S subset C), and
% E_approx=E(empty)*product(K(C), retained nonempty C).
% For noise, K(C)=E(C)-E(empty)-sum(K(S), proper nonempty S),
% and E_approx=E(empty)+sum(K(C)). Masks remain explicit;
% unresolved denominators are never replaced by epsilon.
%
% See also cce_clusters, clexpn.
% Reference: https://doi.org/10.1103/PhysRevB.78.085315
%
% ilya.kuprov@weizmann.ac.il

function reduced=cce_reduce(raw,clusters,parameters)

% Validate indexed raw records and the closed cluster family
grumble(raw,clusters,parameters);

% Native ensemble reduction retains experiment-owned output axes
if isfield(parameters,'cce')

    options=parameters.cce;
    count=numel(raw);
    weights=zeros(1,count);
    for n=1:count
        weights(n)=raw{n}(1).weight;
    end
    selection=struct('method','gcce','order',options.order,...
                     'reducer',options.reducer,'numerics',options.numerics,'expansion',options.expansion);
    selection.numerics.zero_policy='mask';
    if isfield(selection.numerics,'propagator')
        selection.numerics=rmfield(selection.numerics,'propagator');
    end
    additive=isfield(parameters,'noise');
    if additive
        selection.method='noise';
    end

    % Centre every conditional moment with the ensemble means before inclusion-exclusion
    if isfield(parameters,'noise')&&parameters.noise.connected
        for k=1:numel(raw{1})
            mean_time=zeros(size(raw{1}(k).mean_time));
            mean_zero=0;
            for n=1:count
                mean_time=mean_time+weights(n)*raw{n}(k).mean_time;
                mean_zero=mean_zero+weights(n)*raw{n}(k).mean_zero;
            end
            for n=1:count
                raw{n}(k).value=raw{n}(k).moment-raw{n}(k).mean_time*mean_zero-...
                                mean_time*(raw{n}(k).mean_zero-mean_zero);
                raw{n}(k).valid=raw{n}(k).valid&isfinite(mean_time)&isfinite(mean_zero);
            end
        end
    end

    % Locate the exact record on the selected domain, including the empty bath
    bath=sort(options.bath_spins(:));
    sites=[{[]} {clusters.sites}];
    complete=find(cellfun(@(ids)isequal(sort(ids(:)),bath),sites),1);
    interlaced=strcmp(options.bath_state.mode,'interlaced');
    if interlaced&&any([raw{1}.full_bath])
        complete=find([raw{1}.full_bath],1);
    end
    sampling_error=[];

    % The internal-first estimator includes the observed-spin-only reference
    if strcmp(options.averaging,'internal_then_expansion')&&count>1
        records=raw{1};
        for k=1:numel(records)
            records(k).value(:)=0;
            records(k).abs_error_bound(:)=0;
            records(k).valid(:)=true;
            for n=1:count
                if weights(n)==0
                    continue;
                end
                records(k).value=records(k).value+weights(n)*raw{n}(k).value;
                records(k).abs_error_bound=records(k).abs_error_bound+...
                                          weights(n)*raw{n}(k).abs_error_bound;
                records(k).valid=records(k).valid&raw{n}(k).valid;
            end
        end
        sampling_error=NaN(size(records(1).value));
        if ~isempty(complete)&&~(additive&&isfield(parameters,'noise')&&parameters.noise.connected)
            values=[];
            for n=1:count
                values=cat(3,values,raw{n}(complete).value);
            end
            sampling_error=sqrt(sum(abs(values-mean(values,3)).^2,3)/(count*(count-1)));
        end
        raw={records};
        weights=1;
        count=1;
    end

    % Interlaced records already contain parent-conditioned irreducible factors
    parts=cell(1,count);
    for n=1:count
        if interlaced
            parts{n}=interlaced_product(raw{n},clusters,selection);
        else
            parts{n}=cce_reduce(raw{n},clusters,selection);
        end
        if ~isempty(complete)
            exact=raw{n}(complete);
            parts{n}.value=exact.value;
            parts{n}.valid=exact.valid;
            parts{n}.abs_error_bound=exact.abs_error_bound;
            parts{n}.by_order(:,:,end)=exact.value;
            parts{n}.by_order_valid(:,:,end)=exact.valid;
            parts{n}.by_order_error(:,:,end)=exact.abs_error_bound;
            parts{n}.reason(:)={'ok'};
            parts{n}.reason(~exact.valid)={'raw_invalid'};
            parts{n}.by_order_reason(:,:,end)=parts{n}.reason;
        end
    end
    reduced=parts{1};
    reduced.value(:)=0;
    reduced.abs_error_bound(:)=0;
    reduced.valid(:)=true;
    reduced.by_order(:)=0;
    reduced.by_order_error(:)=0;
    reduced.by_order_valid(:)=true;
    for n=1:count
        if weights(n)==0
            continue;
        end
        reduced.value=reduced.value+weights(n)*parts{n}.value;
        reduced.abs_error_bound=reduced.abs_error_bound+weights(n)*parts{n}.abs_error_bound;
        reduced.valid=reduced.valid&parts{n}.valid;
        reduced.by_order=reduced.by_order+weights(n)*parts{n}.by_order;
        reduced.by_order_error=reduced.by_order_error+weights(n)*parts{n}.by_order_error;
        reduced.by_order_valid=reduced.by_order_valid&parts{n}.by_order_valid;
    end

    % Retain the deterministic weighted-summation estimate separately from sampling
    if count>1
        operations=4*count;
        rounding=operations*eps/(1-operations*eps);
        magnitude=zeros(size(reduced.value));
        for n=1:count
            magnitude=magnitude+weights(n)*abs(parts{n}.value);
        end
        reduced.abs_error_bound=reduced.abs_error_bound+rounding*magnitude;
    end

    % A finite exact reference remains valid even above the requested error budget
    reduced.valid=reduced.valid&isfinite(reduced.value);
    if isempty(complete)
        reduced.valid=reduced.valid&reduced.abs_error_bound<=options.numerics.atol+...
                                  options.numerics.rtol*abs(reduced.value);
    end
    reduced.reason=repmat({'ok'},size(reduced.value));
    reduced.reason(~reduced.valid)={'ensemble_invalid'};
    reduced.by_order_reason=repmat({'ok'},size(reduced.by_order));
    reduced.by_order_reason(~reduced.by_order_valid)={'ensemble_invalid'};
    if strcmp(options.bath_state.mode,'interlaced')&&any([raw{1}.full_bath])
        reduced.by_order_reason(:,:,1:end-1)={'not_evaluated'};
    end
    if strcmp(options.numerics.zero_policy,'error')&&any(~reduced.valid,'all')
        error('cce_reduce:unresolved','The native expansion contains unresolved output elements.');
    end
    reduced.samples=parts;
    reduced.weights=weights;
    reduced.complete_record=complete;
    reduced.sampling_standard_error=zeros(size(reduced.value));
    if strcmp(options.bath_state.mode,'sampled')
        reduced.sampling_standard_error(:)=NaN;
    end
    if strcmp(options.bath_state.mode,'sampled')&&count>1&&all(abs(weights-1/count)<=eps)
        values=cellfun(@(part)part.value,parts,'UniformOutput',false);
        values=cat(3,values{:});
        reduced.sampling_standard_error=sqrt(sum(abs(values-mean(values,3)).^2,3)/(count*(count-1)));
    end
    if strcmp(options.bath_state.mode,'interlaced')
        reduced.sampling_standard_error(:)=NaN;
        if count>1
            values=cellfun(@(part)part.value,parts,'UniformOutput',false);
            reduced.sampling_standard_error=std(cat(3,values{:}),0,3)/sqrt(count);
        end
        if any([raw{1}.full_bath])
            reduced.sampling_standard_error(:)=0;
        end
    end
    if ~isempty(sampling_error)
        reduced.sampling_standard_error=sampling_error;
    end
    return

end

% Canonicalise gathered results by stable cluster IDs
[~,baseline_index]=ismember(uint64(0),[raw.cluster_id]);
baseline=raw(baseline_index);
nc=numel(clusters);
[~,indices]=ismember([clusters.id],[raw.cluster_id]);
records=raw(indices);
shape=size(baseline.value);
max_order=max([0 clusters.order]);
by_order=complex(nan([shape max_order+1]));
errors=inf([shape max_order+1]);
valid=false([shape max_order+1]);
reasons=repmat({'ok'},[shape max_order+1]);
reducer='powers';
if isfield(parameters,'reducer')
    reducer=parameters.reducer;
end

% Reduce each prefix independently; powers change with the included family
for o=0:max_order

    % Compile exact signed powers of this order's downward-closed family
    active=find([clusters.order]<=o);
    powers=zeros(1,nc,'int64');
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
    values=complex(nan(shape));
    bounds=inf(shape);
    masks=false(shape);
    why=repmat({'ok'},shape);
    for q=1:numel(values)

        % Gather absolute input uncertainties and validity
        z=complex(zeros(1,nc+1));
        eta=zeros(1,nc+1);
        usable=false(1,nc+1);
        z(1)=baseline.value(q);
        eta(1)=baseline.abs_error_bound(q);
        usable(1)=baseline.valid(q);
        for k=active
            z(k+1)=records(k).value(q);
            eta(k+1)=records(k).abs_error_bound(q);
            usable(k+1)=records(k).valid(q);
        end

        % Empty-bath/order-zero result is the unmodified observed baseline
        if isempty(active)
            values(q)=z(1);
            bounds(q)=eta(1);
            masks(q)=usable(1);
            if ~usable(1)
                why{q}='raw_invalid';
            elseif bounds(q)>parameters.numerics.atol+...
                             parameters.numerics.rtol*abs(values(q))
                masks(q)=false;
                why{q}='error_budget';
            end
            continue
        end

        % Additive connected quantities use signed coefficients, not ratios
        if strcmp(parameters.method,'noise')
            terms=[1-sum(double(powers)) double(powers)];
            needed=find(terms~=0);
            if any(~usable(needed))
                why{q}='raw_invalid';
                continue
            end
            total=complex(0);
            correction=complex(0);
            magnitude=0;
            uncertainty=0;
            for k=needed
                term=terms(k)*z(k);
                increment=term-correction;
                updated=total+increment;
                correction=(updated-total)-increment;
                total=updated;
                magnitude=magnitude+abs(term);
                uncertainty=uncertainty+abs(terms(k))*eta(k);
            end
            operations=8*numel(needed);
            unit=eps/2;
            if operations*unit>=1
                why{q}='arithmetic_range';
                continue
            end
            values(q)=total;
            bounds(q)=uncertainty+operations*unit/(1-operations*unit)*magnitude;
            if ~isfinite(total)||~isfinite(bounds(q))
                why{q}='arithmetic_range';
                continue
            end
            masks(q)=true;
        else

            % Normalisation always requires a resolved nonzero baseline
            if ~usable(1)
                why{q}='raw_invalid';
                continue
            elseif abs(z(1))<=eta(1)
                why{q}='baseline_zero';
                continue
            end

            % Normalise once and propagate the quotient uncertainty
            g=complex(zeros(1,nc));
            ge=zeros(1,nc);
            needed=find(powers~=0);
            if strcmp(reducer,'direct')
                needed=active;
            end
            if any(~usable(needed+1))
                why{q}='raw_invalid';
                continue
            end
            unit=eps/2;
            gamma=6*unit/(1-6*unit);
            for k=needed
                g(k)=z(k+1)/z(1);
                ge(k)=(eta(k+1)+abs(g(k))*eta(1))/(abs(z(1))-eta(1))+gamma*abs(g(k));
            end

            % Use compensated signed log products without phase unwrapping
            if strcmp(reducer,'powers')
                [total,bound,reason]=power_product(z(1),eta(1),g,ge,powers);
                if ~strcmp(reason,'ok')
                    why{q}=reason;
                    continue
                end
            else

                % Evaluate direct irreducible factors in subset order
                factors=complex(ones(1,nc));
                factor_errors=zeros(1,nc);
                reason='ok';
                for k=active
                    subsets=clusters(k).proper_subsets;
                    [denominator,den_error,reason]=power_product(1,0,factors(subsets),...
                                                                factor_errors(subsets),...
                                                                ones(size(subsets),'int64'));
                    if ~strcmp(reason,'ok')
                        break;
                    end
                    if abs(denominator)<=den_error
                        reason='denominator_unresolved';
                        break
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
                    why{q}=reason;
                    continue
                end

            end
            values(q)=total;
            bounds(q)=bound;
            masks(q)=true;

        end

        % Reject estimates exceeding the caller's absolute/relative budget
        if masks(q)&&bounds(q)>parameters.numerics.atol+...
                               parameters.numerics.rtol*abs(values(q))
            masks(q)=false;
            why{q}='error_budget';
        end

    end

    % Preserve invalid masks and reason codes, never plausible finite fillers
    values(~masks)=NaN;
    by_order(:,:,o+1)=values;
    errors(:,:,o+1)=bounds;
    valid(:,:,o+1)=masks;
    reasons(:,:,o+1)=why;

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

% Allocate compensated amplitude, phase, and uncertainty accumulators
value=NaN;
bound=Inf;
reason='ok';
indices=find(powers~=0);
amplitude=log(abs(base));
phase=angle(base);
ac=0;
pc=0;
condition=base_error/(abs(base)-base_error);
exact_zero=false;
scale=abs(amplitude)+abs(phase);
count=0;

% Inspect all denominators before accepting any exact numerator zero
for k=indices

    % Refuse nonfinite intermediate quotients before accepting exact zeros
    if ~isfinite(z(k))||~isfinite(eta(k))
        reason='arithmetic_range';
        return
    end
    if powers(k)<0&&abs(z(k))<=eta(k)
        reason='denominator_unresolved';
        return
    elseif abs(z(k))<=eta(k)
        if z(k)==0&&eta(k)==0&&powers(k)>0
            exact_zero=true;
            continue
        else
            reason='uncertain_zero';
            return
        end
    end
    if abs(powers(k))>flintmax
        reason='arithmetic_range';
        return
    end

    % Integer powers preserve complex phase without logarithm manifold fixes
    p=double(powers(k));
    a=p*log(abs(z(k)));
    b=p*angle(z(k));
    increment=a-ac;
    updated=amplitude+increment;
    ac=(updated-amplitude)-increment;
    amplitude=updated;
    increment=b-pc;
    updated=phase+increment;
    pc=(updated-phase)-increment;
    phase=updated;
    condition=condition+abs(p)*eta(k)/(abs(z(k))-eta(k));
    scale=scale+abs(a)+abs(b);
    count=count+1;

end

% Exact positive zeros need no logarithm or epsilon replacement
if exact_zero
    value=0;
    bound=0;
    return
end

% Convert conditioning and rounding bounds before exponentiating
operations=24*count+8;
unit=eps/2;
if operations*unit>=1||~isfinite(amplitude)||~isfinite(phase)
    reason='arithmetic_range';
    return
end
condition=condition+operations*unit/(1-operations*unit)*(1+scale);
if amplitude>log(realmax)||amplitude<log(realmin*eps)||~isfinite(condition)
    reason='arithmetic_range';
    return
end
value=exp(amplitude)*complex(cos(phase),sin(phase));
if condition==0
    bound=0;
else
    log_error=condition+log(-expm1(-condition));
    if amplitude+log_error>log(realmax)
        reason='arithmetic_range';
        value=NaN;
        return
    end
    bound=exp(amplitude+log_error);
end

end

% Multiply each irreducible factor once, without feeding it into a raw reducer
function reduced=interlaced_product(raw,clusters,selection)
empty=selection;
empty.order=0;
baseline=raw(1);
complete=find([raw.full_bath],1);
if ~isempty(complete)
    baseline=raw(complete);
    baseline.cluster_id=0;
end
reduced=cce_reduce(baseline,clusters([]),empty);
orders=0:selection.order;
shape=[size(raw(1).value) numel(orders)];
reduced.by_order=complex(NaN(shape));
reduced.by_order_error=Inf(shape);
reduced.by_order_valid=false(shape);
reduced.orders=orders;
reduced.by_order_reason=repmat({'not_evaluated'},shape);
if ~isempty(complete)
    reduced.by_order(:,:,end)=baseline.value;
    reduced.by_order_error(:,:,end)=baseline.abs_error_bound;
    reduced.by_order_valid(:,:,end)=baseline.valid;
    reduced.by_order_reason(:,:,end)=reduced.reason;
    return
end
value=baseline.value;
bound=baseline.abs_error_bound;
valid=baseline.valid;
for o=orders
    for k=find([clusters.order]==o)
        factor=raw(k+1);
        bound=abs(value).*factor.factor_abs_error_bound+...
              abs(factor.factor).*bound+bound.*factor.factor_abs_error_bound;
        value=value.*factor.factor;
        valid=valid&factor.factor_valid;
    end
    reduced.by_order(:,:,o+1)=value;
    reduced.by_order_error(:,:,o+1)=bound;
    reduced.by_order_valid(:,:,o+1)=valid;
    reasons=repmat({'ok'},size(valid));
    reasons(~valid)={'raw_invalid'};
    reduced.by_order_reason(:,:,o+1)=reasons;
end
reduced.value=value;
reduced.abs_error_bound=bound;
reduced.valid=valid;
end

% Validate raw gather indexing, scalar budgets and the subset order
function grumble(raw,clusters,parameters)
if ~isstruct(parameters)||~isscalar(parameters)
    error('cce_reduce:parameters','method and numerical budgets required.');
end
if isfield(parameters,'cce')
    options=parameters.cce;
    if ~iscell(raw)||isempty(raw)||~isvector(raw)||~isstruct(clusters)||...
       ~isstruct(options)||~isscalar(options)||...
       ~all(isfield(options,{'order','reducer','numerics','expansion','averaging','bath_spins','bath_state'}))
        error('cce_reduce:inputs','native raw cells and resolved cluster controls required.');
    end
    count=numel(raw);
    weights=zeros(1,count);
    for n=1:count
        if ~isstruct(raw{n})||isempty(raw{n})||~isfield(raw{n},'weight')||...
           ~isnumeric(raw{n}(1).weight)||~isreal(raw{n}(1).weight)||~isscalar(raw{n}(1).weight)
            error('cce_reduce:weights','Configuration weights must be normalized and nonnegative.');
        end
        weights(n)=raw{n}(1).weight;
    end
    if any(~isfinite(weights))||any(weights<0)||abs(sum(weights)-1)>count*eps
        error('cce_reduce:weights','Configuration weights must be normalized and nonnegative.');
    end
    if strcmp(options.bath_state.mode,'interlaced')
        if isfield(parameters,'noise')
            error('cce_reduce:interlaced','Interlacing requires multiplicative signals.');
        end
        if ~strcmp(options.averaging,'after_expansion')
            error('cce_reduce:interlaced','Interlaced factors require after_expansion averaging.');
        end
    end
    selection=struct('method','gcce','reducer',options.reducer,'numerics',options.numerics);
    for n=1:count
        records=raw{n};
        family=clusters;
        if strcmp(options.bath_state.mode,'interlaced')&&any([records.full_bath])
            records=records(find([records.full_bath],1));
            records.cluster_id=0;
            family=clusters([]);
        end
        grumble(records,family,selection);
    end
    return
end
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


