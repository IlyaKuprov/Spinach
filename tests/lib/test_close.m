% Adds a numerical regression check with tolerances and explanation. Syntax:
%
%       result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)
%
% Parameters:
%
%     result     - test result structure
%
%     label      - check label
%
%     observed   - value produced by Spinach
%
%     reference  - independently known right answer
%
%     abs_tol    - absolute tolerance
%
%     rel_tol    - relative tolerance
%
%     why        - explanation of the right answer
%
% Outputs:
%
%     result     - updated test result structure
%
% A failed check is recorded in the messages and failures fields of the
% result structure rather than thrown, so that the checks that follow it
% in the calling test are still evaluated and the earlier passes are not
% lost; run_tests inspects the failures field to decide the status. The
% conditions that prevent the comparison from being made at all, namely
% a size mismatch and non-finite data, are recorded in the same way.
%
% ilya.kuprov@weizmann.ac.il

function result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)

% Convert sparse arrays for norm evaluation
observed=full(observed);
reference=full(reference);

% Convert data to double vectors for numerical comparison
observed_vec=double(observed(:));
reference_vec=double(reference(:));

% Identify a condition that prevents the comparison
if ~isequal(size(observed),size(reference))
    detail=[label ' -- size mismatch'];
elseif any(~isfinite(observed_vec))
    detail=[label ' -- observed value contains NaN or Inf'];
elseif any(~isfinite(reference_vec))
    detail=[label ' -- reference value contains NaN or Inf'];
else
    detail='';
end

% Record an impossible comparison and return
if ~isempty(detail)
    result.messages{end+1}=['FAIL: ' detail];
    result.failures{end+1}=detail;
    return
end

% Compute a scaled Frobenius/vector norm error
error_norm=norm(observed_vec-reference_vec,2);
ref_norm=max([1 norm(reference_vec,2)]);
limit=abs_tol+rel_tol*ref_norm;

% Identify a non-finite comparison or an error above the tolerance
if (~isfinite(error_norm))||(~isfinite(limit))
    detail=[label ' -- comparison produced a non-finite scalar'];
elseif error_norm>limit
    detail=[label ', error=' num2str(error_norm,'%.3e') ...
            ', limit=' num2str(limit,'%.3e') ' -- ' why];
end

% Record a failed comparison and return
if ~isempty(detail)
    result.messages{end+1}=['FAIL: ' detail];
    result.failures{end+1}=detail;
    return
end

% Record the pass message
result.messages{end+1}=['PASS: ' label ', error=' ...
                        num2str(error_norm,'%.3e') ', tolerance=' ...
                        num2str(limit,'%.3e') ' -- ' why];

end

