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

% Check consistency
grumble(result,label,observed,reference,abs_tol,rel_tol,why);

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

% Compare against the tolerance where the comparison is possible
if isempty(detail)
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
end

% Record the outcome of the check
if isempty(detail)
    result.messages{end+1}=['PASS: ' label ', error=' ...
                            num2str(error_norm,'%.3e') ', tolerance=' ...
                            num2str(limit,'%.3e') ' -- ' why];
else
    result.messages{end+1}=['FAIL: ' detail];
    result.failures{end+1}=detail;
end

% Retain the record for the run_tests catch path
test_record(result);

end

% Consistency enforcement
function grumble(result,label,observed,reference,abs_tol,rel_tol,why)
if (~isstruct(result))||(~isscalar(result))||...
   (~isfield(result,'messages'))||(~isfield(result,'failures'))
    error('result must be a scalar test result structure.');
end
if (~ischar(label))||isempty(label)||(~isrow(label))
    error('label must be a non-empty character string.');
end
if (~islogical(observed))&&(~isnumeric(observed))
    error('observed must be logical or numeric.');
end
if (~islogical(reference))&&(~isnumeric(reference))
    error('reference must be logical or numeric.');
end
if (~isnumeric(abs_tol))||(~isscalar(abs_tol))||(~isreal(abs_tol))||(abs_tol<0)
    error('abs_tol must be a non-negative real scalar.');
end
if (~isnumeric(rel_tol))||(~isscalar(rel_tol))||(~isreal(rel_tol))||(rel_tol<0)
    error('rel_tol must be a non-negative real scalar.');
end
if (~ischar(why))||isempty(why)||(~isrow(why))
    error('why must be a non-empty character string.');
end
end


