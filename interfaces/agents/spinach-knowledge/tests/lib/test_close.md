# tests/lib/test_close.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_close.m`
- Signature: `result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)`
- Total lines: 71

## Purpose

Adds a numerical regression check with tolerances and explanation. Syntax: result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Convert sparse arrays for norm evaluation; implemented by `observed=full(observed)`.
- Lines 33-34: Check dimensions first; implemented by `if ~isequal(size(observed),size(reference))`.
- Lines 38-39: Convert data to double vectors for numerical comparison; implemented by `observed_vec=double(observed(:))`.
- Lines 42-43: Reject non-finite values before norm evaluation; implemented by `if any(~isfinite(observed_vec))`.
- Lines 50-51: Compute a scaled Frobenius/vector norm error; implemented by `error_norm=norm(observed_vec-reference_vec,2)`.
- Lines 55-56: Reject non-finite comparison scalars; implemented by `if (~isfinite(error_norm))||(~isfinite(limit))`.
- Lines 60-61: Fail with the numerical details; implemented by `if error_norm>limit`.
- Lines 66-69: Record the pass message; implemented by `result.messages{end+1}=['PASS: ' label ', error=' num2str(error_norm,'%.3e') ', tolerance=' num2str(limit,'%.3e') ' -- ' why]`.

### Control flow inferred from the code

- Line 34: conditional branch on `~isequal(size(observed),size(reference))`.
- Line 43: conditional branch on `any(~isfinite(observed_vec))`.
- Line 46: conditional branch on `any(~isfinite(reference_vec))`.
- Line 56: conditional branch on `(~isfinite(error_norm))||(~isfinite(limit))`.
- Line 61: conditional branch on `error_norm>limit`.

### Key state/data transformations

- Lines 30: computes `observed` using `observed=full(observed)`.
- Lines 31: computes `reference` using `reference=full(reference)`.
- Lines 39: computes `observed_vec` using `observed_vec=double(observed(:))`.
- Lines 40: computes `reference_vec` using `reference_vec=double(reference(:))`.
- Lines 51: computes `error_norm` using `error_norm=norm(observed_vec-reference_vec,2)`.
- Lines 52: computes `ref_norm` using `ref_norm=max([1 norm(reference_vec,2)])`.
- Lines 53: computes `limit` using `limit=abs_tol+rel_tol*ref_norm`.
- Lines 62-63: computes `error(['FAILED: ' label ', error` using `error(['FAILED: ' label ', error=' num2str(error_norm,'%.3e') ', limit=' num2str(limit,'%.3e') ' -- ' why])`.
- Lines 63: computes `', limit` using `', limit=' num2str(limit,'%.3e') ' -- ' why])`.
- Lines 67-69: computes `result.messages{end+1}` using `result.messages{end+1}=['PASS: ' label ', error=' num2str(error_norm,'%.3e') ', tolerance=' num2str(limit,'%.3e') ' -- ' why]`.
- Lines 68-69: computes `num2str(error_norm,'%.3e') ', tolerance` using `num2str(error_norm,'%.3e') ', tolerance=' num2str(limit,'%.3e') ' -- ' why]`.

## Parameters / inputs

- result -test result structure
- label -check label
- observed -value produced by Spinach
- reference -independently known right answer
- abs_tol -absolute tolerance
- rel_tol -relative tolerance
- why -explanation of the right answer

## Outputs

- result -updated test result structure

## Implementation structure

- Adds a numerical regression check with tolerances and explanation. Syntax:
- result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)
- result -test result structure
- label -check label
- observed -value produced by Spinach
- reference -independently known right answer
- abs_tol -absolute tolerance
- rel_tol -relative tolerance
- why -explanation of the right answer
- result -updated test result structure
- Convert sparse arrays for norm evaluation
- Check dimensions first

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isequal()`, `double()`, `observed()`, `reference()`, `any()`, `num2str()`.
