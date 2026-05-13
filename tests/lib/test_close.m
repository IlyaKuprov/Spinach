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
% ilya.kuprov@weizmann.ac.il

function result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)

% Convert sparse arrays for norm evaluation
observed=full(observed);
reference=full(reference);

% Check dimensions first
if ~isequal(size(observed),size(reference))
    error(['FAILED: ' label ' -- size mismatch.']);
end

% Convert data to double vectors for numerical comparison
observed_vec=double(observed(:));
reference_vec=double(reference(:));

% Reject non-finite values before norm evaluation
if any(~isfinite(observed_vec))
    error(['FAILED: ' label ' -- observed value contains NaN or Inf.']);
end
if any(~isfinite(reference_vec))
    error(['FAILED: ' label ' -- reference value contains NaN or Inf.']);
end

% Compute a scaled Frobenius/vector norm error
error_norm=norm(observed_vec-reference_vec,2);
ref_norm=max([1 norm(reference_vec,2)]);
limit=abs_tol+rel_tol*ref_norm;

% Reject non-finite comparison scalars
if (~isfinite(error_norm))||(~isfinite(limit))
    error(['FAILED: ' label ' -- comparison produced a non-finite scalar.']);
end

% Fail with the numerical details
if error_norm>limit
    error(['FAILED: ' label ', error=' num2str(error_norm,'%.3e') ...
           ', limit=' num2str(limit,'%.3e') ' -- ' why]);
end

% Record the pass message
result.messages{end+1}=['PASS: ' label ', error=' ...
                        num2str(error_norm,'%.3e') ', tolerance=' ...
                        num2str(limit,'%.3e') ' -- ' why];

end
