% Matrix-free core for a polyadic, with forward and adjoint block actions.
% Syntax:
%
%       core=matfree(dims,forward,adjoint,real_flag)
%
% Parameters:
%
%    dims       - [rows columns], both at least two
%    forward    - handle mapping columns-by-n blocks to rows-by-n blocks
%    adjoint    - handle mapping rows-by-n blocks to columns-by-n blocks
%    real_flag  - logical scalar, true if the represented matrix is real
%
% Outputs:
%
%    core       - implicit matrix supporting polyadic exponential actions
%
% Actions must be linear, finite, mutually adjoint, and preserve the device
% of their input. Their internal data are the caller's responsibility.
% Scalar multiplication scales the operator, not a one-column state.
% Explicit materialisation is deliberately unavailable.
%
% ilya.kuprov@weizmann.ac.il

classdef (InferiorClasses={?gpuArray}) matfree

    % Immutable operator description
    properties (SetAccess=private)
        dims
        forward
        adjoint
        real_flag
        coeff=1;
    end

    methods

        % Store the actions and their matrix dimensions
        function core=matfree(dims,forward,adjoint,real_flag)
            grumble(dims,forward,adjoint,real_flag);
            core.dims=dims; core.forward=forward;
            core.adjoint=adjoint; core.real_flag=real_flag;
        end

        % Matrix compatibility for polyadic validation
        function answer=isnumeric(core) %#ok<MANU>
            answer=true;
        end

        % Matrix compatibility for Kronecker contractions
        function answer=ismatrix(core) %#ok<MANU>
            answer=true;
        end

        % Matrix element count for scalar dispatch
        function answer=numel(core)
            answer=prod(core.dims);
        end

        % An implicit action is never declared to be an identity
        function answer=iseye(core) %#ok<MANU>
            answer=false;
        end

        % Structural nonzero marker, not an entry count
        function answer=nnz(core)
            answer=double(core.coeff~=0);
        end

        % Finiteness of the stored coefficient, assuming valid actions
        function answer=allfinite(core)
            answer=isfinite(core.coeff);
        end

        % Reality of the represented operator
        function answer=isreal(core)
            answer=core.real_flag&&isreal(core.coeff);
        end

        % Actions receive GPU blocks without transforming their closures
        function core=gpuArray(core)
        end

        % Explicit matrices cannot be recovered from opaque actions
        function answer=full(core) %#ok<STOUT,MANU>
            error('matrix-free cores cannot be materialised.');
        end

    end

end

% Consistency enforcement
function grumble(dims,forward,adjoint,real_flag)
validateattributes(dims,{'double'},{'row','numel',2,'integer','finite','>=',2});
if ~isa(forward,'function_handle')||~isa(adjoint,'function_handle')
    error('forward and adjoint must be function handles.');
end
if ~islogical(real_flag)||~isscalar(real_flag)
    error('real_flag must be a logical scalar.');
end
end


