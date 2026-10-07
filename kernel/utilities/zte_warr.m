% Leakage-based time estimate for zero track elimination. Syntax:
%
%          duration=zte_warr(spin_system,L,rho,projector)
%
% Parameters:
%
%    spin_system - Spinach structure with the positive absolute
%                  state 2-norm tolerance tols.zte_warr
%
%    L           - fixed square Liouvillian, in radians per second
%
%    rho         - initial state column vector
%
%    projector   - coordinate embedding returned by zte, or the
%                  scalar 1 when no reduction is performed
%
% Output:
%
%    duration    - estimated time in seconds before the absolute
%                  state 2-norm error reaches tols.zte_warr
%
% Notes: initial discarded norm d consumes the error tolerance. The
%        leakage rate r is the discarded norm of L acting on the
%        retained initial state. The local estimate (tolerance-d)/r
%        is capped at 1/cheap_norm(L), the ZTE exploration step.
%        The cap prevents extrapolation beyond this short probe
%        scale when initial leakage vanishes or cancels; it is NOT
%        an error bound. Later leakage and amplification are not
%        controlled, even for unitary dynamics. Zero initial rate
%        gives the exploration step, not infinite validity.
%        Initial discard at or above tolerance gives zero; no
%        reduction or a zero input state gives Inf. A zero generator
%        with nonzero state gives NaN (no informative leakage time).
%        Only the supplied vector and fixed generator are covered,
%        not spectra, numerical propagation error, or roundoff.
%        L must be finite; only its action and cheap norm are checked.
%        Apart from cheap_norm (CPU 1-norm, GPU infinity-norm),
%        only one matrix-vector product and vector operations occur.
%
% ilya.kuprov@weizmann.ac.il

function duration=zte_warr(spin_system,L,rho,projector)

% Validate the input
grumble(spin_system,L,rho,projector);

% Identify the retained coordinates without expanding a scalar identity
if isequal(projector,1)
    retained=true(size(rho));
else
    retained=any(projector,2);
end

% Account for initial discard and exact unchanged-state cases
discarded=norm(rho(~retained)); duration=NaN;
if discarded>=spin_system.tols.zte_warr
    duration=0;
elseif all(retained)||~any(rho)
    duration=Inf;
else

    % Evaluate the initial leakage using only the retained state
    rho(~retained)=0; action=L*rho;
    rate=norm(action(~retained)); scale=cheap_norm(L);
    if any(~isfinite(action))||~isfinite(scale)
        error('non-finite ZTE leakage probe; L must be finite.');
    end

    % Cap the local extrapolation at the ZTE exploration step
    if scale>0
        duration=min((spin_system.tols.zte_warr-discarded)/rate,1/scale);
    end
end

% Report the estimate without claiming an error guarantee
report(spin_system,['ZTE warranty estimate: absolute 2-norm tolerance ' ...
                    num2str(spin_system.tols.zte_warr,17) ', estimated time ' ...
                    num2str(duration,17) ' seconds (supplied vector, fixed generator; not a guarantee).']);

end

% Input validation function
function grumble(spin_system,L,rho,projector)
if (~isnumeric(L))||(~ismatrix(L))||(size(L,1)~=size(L,2))
    error('L must be a square numeric matrix.');
end
if (~isnumeric(rho))||(size(rho,2)~=1)||(size(rho,1)~=size(L,1))||...
   any(~isfinite(nonzeros(rho)))
    error('rho must be a finite column vector matching L.');
end
if (~isnumeric(projector))||(~ismatrix(projector))||...
   (~isequal(projector,1)&&(size(projector,1)~=size(L,1)))
    error('projector must be the coordinate embedding returned by zte.');
end
if (~isnumeric(spin_system.tols.zte_warr))||(~isreal(spin_system.tols.zte_warr))||...
   (~isscalar(spin_system.tols.zte_warr))||(~isfinite(spin_system.tols.zte_warr))||...
   (spin_system.tols.zte_warr<=0)
    error('zte_warr must be a finite positive real scalar.');
end
end


