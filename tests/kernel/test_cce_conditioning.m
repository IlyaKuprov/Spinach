% Controlled CCE zeros, uncertainty, cancellation and phase tests. Syntax:
%
%                    result=test_cce_conditioning()
%
% Outputs:
%    result - regression result with exact mask and reason checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_conditioning()

% Set up a path pair expansion with a negative central-vertex exponent
result=new_test_result('kernel/cce_conditioning','CCE zero conditioning',...
                       'Expose unresolved ratios rather than regularise or lose phase.');
model=struct('physical_ids',1:4,'central_ids',1,'bath_ids',2:4,...
             'multiplicities',[2 2 2 2]);
adjacency=logical([0 1 0;1 0 1;0 1 0]);
graph=struct('bath_ids',2:4,'adjacency',sparse(adjacency),'scores',double(adjacency));
parameters=struct('method','cce','order',2,'expansion','standard',...
                  'numerics',struct('atol',1e-7,'rtol',1e-7));
clusters=cce_clusters(model,graph,parameters);
raw=repmat(struct('cluster_id',uint64(0),'value',ones(6,1),...
                  'abs_error_bound',zeros(6,1),'valid',true(6,1)),numel(clusters)+1,1);
for k=1:numel(clusters)
    raw(k+1).cluster_id=clusters(k).id;
end

% Plant a baseline zero, unresolved denominator, uncertain zero and exact zero
raw(1).value(2)=0;
raw(3).value(3)=1e-10; raw(3).abs_error_bound(3)=1e-9;
raw(5).value(4)=0; raw(5).abs_error_bound(4)=1e-9;
raw(5).value(5)=0;
raw(5).value(6)=0; raw(3).value(6)=0;
reduced=cce_reduce(raw,clusters,parameters);
expected={'ok';'baseline_zero';'denominator_unresolved';'uncertain_zero';'ok';'denominator_unresolved'};
result=test_true(result,'Exact conditioning masks and reason codes',...
                 isequal(reduced.valid,[true;false;false;false;true;false])&&...
                 isequal(reduced.reason,expected)&&reduced.value(5)==0&&...
                 all(isnan(reduced.value(~reduced.valid))),...
                 'six controlled cases distinguish baseline, numerator and denominator uncertainty');

% An unresolved denominator wins over an exact positive zero elsewhere
parameters.numerics.zero_policy='error'; rejected=false;
try
    cce_reduce(raw,clusters,parameters);
catch exception
    rejected=strcmp(exception.identifier,'cce_reduce:invalid');
end
result=test_true(result,'Explicit fail-fast zero policy',rejected,...
                 'error policy refuses unresolved final values instead of silently masking them');
parameters.numerics.zero_policy='mask';

% Preserve a signed phase crossing pi with no unwrapping through zeros
for k=1:numel(raw)
    raw(k).value=ones(6,1); raw(k).abs_error_bound=zeros(6,1);
end
phase=[0;pi-0.01;pi+0.01;2*pi-0.01;2*pi+0.01;-pi-0.01];
raw(3).value=exp(1i*phase); raw(5).value=0.8*exp(2i*phase);
raw(6).value=0.7*exp(-0.4i*phase);
reduced=cce_reduce(raw,clusters,parameters); expected=0.56*exp(0.6i*phase);
result=test_close(result,'Signed phases across logarithm cuts',reduced.value,expected,3e-14,3e-14,...
                  'integer exponents preserve actual complex phase without branch manipulation');

% Negative denominator uncertainty produces an observable error budget failure
raw(3).abs_error_bound(:)=1e-4; reduced=cce_reduce(raw,clusters,parameters);
result=test_true(result,'Caller error-budget refusal',all(~reduced.valid)&&...
                 all(strcmp(reduced.reason,'error_budget')),...
                 'resolved but ill-conditioned products are not reported as accurate');

% Compensated logarithms avoid overflow in intermediate products
for k=1:numel(raw)
    raw(k).value=ones(6,1); raw(k).abs_error_bound=zeros(6,1);
end
raw(3).value(:)=1e300; raw(5).value(:)=1e200; raw(6).value(:)=1e100;
wide=cce_reduce(raw,clusters,parameters);
result=test_close(result,'Resolved logarithmic dynamic range',wide.value,ones(6,1),3e-13,0,...
                  '1e200*1e100/1e300 is reduced before overflowing any intermediate product');

% Final overflow and uncertain baseline are explicitly distinguishable
raw(5).value(:)=1e300; raw(6).value(:)=1e300; raw(3).value(:)=1;
overflow=cce_reduce(raw,clusters,parameters);
result=test_true(result,'Overflow reason',all(~overflow.valid)&&...
                 all(strcmp(overflow.reason,'arithmetic_range')),...
                 'unrepresentable final magnitude is not returned as an infinity');
raw(1).value(:)=1e-10; raw(1).abs_error_bound(:)=1e-9;
baseline=cce_reduce(raw,clusters,parameters);
result=test_true(result,'Uncertain nonzero baseline',all(~baseline.valid)&&...
                 all(strcmp(baseline.reason,'baseline_zero')),...
                 'baseline resolution uses its uncertainty, not equality with floating-point zero');

% Cancellation may remove algebraically irrelevant raw zeros in powers mode
parameters.clusters={1:3}; parameters.order=3;
clusters=cce_clusters(model,graph,parameters);
raw=repmat(struct('cluster_id',uint64(0),'value',1,'abs_error_bound',0,'valid',true),...
           numel(clusters)+1,1);
for k=1:numel(clusters)
    raw(k+1).cluster_id=clusters(k).id; raw(k+1).value=0;
end
raw(end).value=0.2+0.3i; reduced=cce_reduce(raw,clusters,parameters);
result=test_close(result,'Exact signed-power cancellation',reduced.value,raw(end).value,8*eps,0,...
                  'full family uses only its maximal raw factor, not zero-power denominators');
parameters.reducer='direct'; direct=cce_reduce(raw,clusters,parameters);
result=test_true(result,'Direct singularity is reported honestly',~direct.valid&&...
                 strcmp(direct.reason{1},'denominator_unresolved'),...
                 'direct ratios cannot divide exact zero irreducibles even when powers cancel them');

end
