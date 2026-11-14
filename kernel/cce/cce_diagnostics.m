% Audit CCE numerical coverage and density physicality without repair.
% Syntax:
%
%       diagnostics=cce_diagnostics(model,raw,result,parameters)
%
% raw is a cell of cluster-record arrays. result supplies value, valid,
% abs_error_bound and optionally density/sampling_standard_error. Bounds
% are numerical estimates, not truncation confidence intervals. Invalid
% reasons, time coverage and density tests are reported independently;
% no normalization, positivity projection or sample censoring is done.
%
% Inputs:
%   result - CCE value, valid mask, abs_error_bound and retained metadata.
%   parameters - numerical tolerances and observable declarations.
%
% Outputs:
%   diagnostics - numerical/conditioning, cluster, sample and density
%   diagnostic fields.
%
% ilya.kuprov@weizmann.ac.il

function diagnostics=cce_diagnostics(model,raw,result,parameters)

% Check shapes before combining numerical validity and time coverage
grumble(model,raw,result,parameters);
value=result.value; valid=result.valid&isfinite(value)&...
      isfinite(result.abs_error_bound)&result.abs_error_bound>=0;
tolerance=parameters.numerics.atol+parameters.numerics.rtol*abs(value);
nt=size(value,1); diagnostics=struct();
diagnostics.valid=valid;
diagnostics.valid_time_mask=all(valid,2);
diagnostics.valid_time_coverage=nnz(diagnostics.valid_time_mask)/nt;
diagnostics.numerical_error=result.abs_error_bound;
diagnostics.numerical_budget_met=valid&result.abs_error_bound<=tolerance;
diagnostics.invalid_reasons=cell(size(value));
for k=1:numel(value)
    reasons=cell(1,4); count=0;
    if ~result.valid(k), count=count+1; reasons{count}='reduction_invalid'; end
    if ~isfinite(value(k)), count=count+1; reasons{count}='nonfinite_value'; end
    if ~isfinite(result.abs_error_bound(k))||result.abs_error_bound(k)<0
        count=count+1; reasons{count}='invalid_error_bound';
    elseif result.abs_error_bound(k)>tolerance(k)
        count=count+1; reasons{count}='numerical_budget';
    end
    diagnostics.invalid_reasons{k}=strjoin(reasons(1:count),',');
end

% Keep raw propagation, reduction and sampling errors distinguishable
diagnostics.maximum_eigen_residual=0;
diagnostics.maximum_orthogonality_residual=0;
diagnostics.raw_invalid=0; diagnostics.raw_records=0;
for s=1:numel(raw)
    for k=1:numel(raw{s})
        record=raw{s}(k); diagnostics.raw_records=diagnostics.raw_records+1;
        diagnostics.raw_invalid=diagnostics.raw_invalid+nnz(~record.valid);
        if isfield(record,'eigen_residual')
            diagnostics.maximum_eigen_residual=max(diagnostics.maximum_eigen_residual,record.eigen_residual);
        end
        if isfield(record,'orthogonality_residual')
            diagnostics.maximum_orthogonality_residual=max(diagnostics.maximum_orthogonality_residual,record.orthogonality_residual);
        end
    end
end
diagnostics.sampling_standard_error=[];
if isfield(result,'sampling_standard_error')
    diagnostics.sampling_standard_error=result.sampling_standard_error;
end
diagnostics.truncation_error=[]; diagnostics.structural_error=[];
diagnostics.sampling_error=[]; diagnostics.time_mesh_error=[];
diagnostics.complete_coverage=all(valid,'all');
diagnostics.physicality_checked=false; diagnostics.physicality_pass=true;

% Inspect density matrices as returned, including all invalid time points
if isfield(result,'density')&&~isempty(result.density)
    density=result.density; dimension=size(density,1);
    diagnostics.physicality_checked=true;
    diagnostics.trace_error=nan(nt,1); diagnostics.hermiticity_error=nan(nt,1);
    diagnostics.minimum_eigenvalue=nan(nt,1); diagnostics.physicality_valid=false(nt,1);
    for k=1:nt
        rho=density(:,:,k);
        if ~all(isfinite(rho),'all'), continue; end
        bound=norm(result.abs_error_bound(k,:));
        threshold=parameters.numerics.atol+parameters.numerics.rtol*norm(rho,'fro');
        diagnostics.trace_error(k)=abs(trace(rho)-1);
        diagnostics.hermiticity_error(k)=norm(rho-rho','fro');
        diagnostics.minimum_eigenvalue(k)=min(real(eig((rho+rho')/2)));
        diagnostics.physicality_valid(k)=diagnostics.valid_time_mask(k)&&...
            diagnostics.trace_error(k)<=sqrt(dimension)*bound+threshold&&...
            diagnostics.hermiticity_error(k)<=2*bound+threshold&&...
            diagnostics.minimum_eigenvalue(k)>=-bound-threshold;
    end
    diagnostics.physicality_pass=all(diagnostics.physicality_valid);
end

% Never promote incomplete coverage or unphysical matrices to success
diagnostics.acceptable=diagnostics.complete_coverage&&...
                       all(diagnostics.numerical_budget_met,'all')&&...
                       diagnostics.physicality_pass;
diagnostics.error_bound_kind='numerical estimate, not CCE truncation bound';

end

% Validate matching result arrays and the explicitly requested time mesh
function grumble(model,raw,result,parameters)
if ~isstruct(model)||~iscell(raw)||~isstruct(result)||...
   ~all(isfield(result,{'value','valid','abs_error_bound'}))||...
   ~isstruct(parameters)||~all(isfield(parameters,{'times','numerics'}))||...
   ~all(isfield(parameters.numerics,{'atol','rtol'}))||...
   isempty(result.value)||~isnumeric(result.value)||...
   ~islogical(result.valid)||~isequal(size(result.value),size(result.valid),size(result.abs_error_bound))||...
   size(result.value,1)~=numel(parameters.times)
    error('cce_diagnostics:input','Result arrays must match the requested time mesh.');
end
for field={'atol','rtol'}
    value=parameters.numerics.(field{1});
    if ~isnumeric(value)||~isscalar(value)||~isreal(value)||~isfinite(value)||value<0
        error('cce_diagnostics:input','Numerical tolerances must be finite nonnegative scalars.');
    end
end
if parameters.numerics.atol+parameters.numerics.rtol==0
    error('cce_diagnostics:input','At least one numerical tolerance must be positive.');
end
if isfield(result,'density')&&~isempty(result.density)&&...
   (size(result.density,1)~=size(result.density,2)||size(result.density,3)~=numel(parameters.times))
    error('cce_diagnostics:input','Density matrices must be square on every requested time.');
end
end

