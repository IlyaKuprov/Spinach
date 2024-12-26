% Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri-
% cal harmonic in the expansion of the product of Y(L1,M1) and Y(L2,M2)
% spherical harmonics. In the more general sense, the coefficient refers
% to the expansion coefficient of |L,M> angular momentum or spin state
% in the product basis of |L1,M1>|L2,M2> states. Syntax:
%
%                   cg=clebsch_gordan(L,M,L1,M1,L2,M2)
%
% Parameters:
%
%     L,M,L1,M1,L2,M2  - integer or half-integer indices of 
%                        the angular momentum or spin states
%
% Outputs:
%
%     cg               - floating-point (double precision) 
%                        Clebsch-Gordan coefficient
%                                
% Note: only some combinations of L,M,L1,M1,L2,M2 are allowed by the pro-
%       perties of spherical harmonics and spin states. If inadmissible
%       indices are supplied, zero is returned.
%
% Note: CG coefficient calculation in double-precision arithmetic is not
%       a trivial matter for high ranks. This function produces machine
%       precision answers up to about L=1e4. A faster implementation for
%       low ranks is available in cg_fast.m function.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=clebsch_gordan.m>

function cg=clebsch_gordan(L,M,L1,M1,L2,M2)

% Check consistency
grumble(L,M,L1,M1,L2,M2);

% Match the notation to Varshalovich, Section 8.2.1
cg=0; c=L; gam=M; a=L1; alp=M1; b=L2; bet=M2;

% Run zero tests (Stage I)
prefactor_is_nonzero=(a+alp>=0)&&(a-alp>=0)&&...
                     (b+bet>=0)&&(b-bet>=0)&&...
                     (c+gam>=0)&&(c-gam>=0)&&...
                     (gam==alp+bet);
    
% Proceed if appropriate
if prefactor_is_nonzero
    
    % Run zero tests (Stage II)
    delta_is_nonzero=(a+b-c>=0)&&(a-b+c>=0)&&...
                     (-a+b+c>=0)&&(a+b+c+1>=0);
                     
    % Proceed if appropriate
    if delta_is_nonzero
            
        % Run zero tests (Stage III)
        lower_sum_limit=max([alp-a, b+gam-a, 0]);
        upper_sum_limit=min([c+b+alp, c+b-a, c+gam]);
        sum_is_nonzero=(upper_sum_limit>=lower_sum_limit);
    
        % Proceed if appropriate
        if sum_is_nonzero
            
            % Compile a look-up table of gamma functions
            gamma(a+b+c+2)=java.math.BigInteger.ONE; gamma(:)=gamma(a+b+c+2);
            for k=3:(a+b+c+2)
                gamma(k)=gamma(k-1).multiply(java.math.BigInteger.valueOf(k-1));
            end
            
            % Compute prefactor numerator
            numer=java.math.BigInteger.valueOf(2*c+1);
            numer=numer.multiply(gamma(a-b+c+1));
            numer=numer.multiply(gamma(-a+b+c+1));
            numer=numer.multiply(gamma(c+gam+1));
            numer=numer.multiply(gamma(c-gam+1));
            numer=numer.multiply(gamma(a+b-c+1));
            
            % Compute prefactor denominator
            denom=gamma(a+b+c+2);
            denom=denom.multiply(gamma(a+alp+1));
            denom=denom.multiply(gamma(a-alp+1));
            denom=denom.multiply(gamma(b+bet+1));
            denom=denom.multiply(gamma(b-bet+1));
            
            % Perform floating-point division
            accur=length(denom.toString)-length(numer.toString)+64;
            numer=java.math.BigDecimal(numer); denom=java.math.BigDecimal(denom);
            prefactor=numer.divide(denom,accur,java.math.RoundingMode.HALF_UP);
            
            % Compute the sum
            z_sum=java.math.BigDecimal.ZERO;
            for z=lower_sum_limit:upper_sum_limit
                
                % Compute term numerator
                numer=java.math.BigInteger.valueOf((-1)^(b+bet+z));
                numer=numer.multiply(gamma(c+b+alp-z+1));
                numer=numer.multiply(gamma(a-alp+z+1));
                
                % Compute term denominator
                denom=gamma(z+1);
                denom=denom.multiply(gamma(c-a+b-z+1));
                denom=denom.multiply(gamma(c+gam-z+1));
                denom=denom.multiply(gamma(a-b-gam+z+1));
                
                % Perform floating-point division
                numer=java.math.BigDecimal(numer); denom=java.math.BigDecimal(denom);
                sum_term=numer.divide(denom,accur,java.math.RoundingMode.HALF_UP);
                
                % Add the term to the total
                z_sum=z_sum.add(sum_term);
                
            end
            
            % Return to double precision
            cg_sq=z_sum.multiply(z_sum).multiply(prefactor);
            cg=z_sum.signum*sqrt(double(cg_sq));
            
        end
            
    end
        
end

end

% Consistency enforcement
function grumble(L,M,L1,M1,L2,M2)
if (~isnumeric(L))||(~isnumeric(L1))||(~isnumeric(L2))||...
   (~isnumeric(M))||(~isnumeric(M1))||(~isnumeric(M2))
    error('all arguments must be numeric.');
end
if (~isreal(L))||(~isreal(L1))||(~isreal(L2))||...
   (~isreal(M))||(~isreal(M1))||(~isreal(M2))
    error('all arguments must be real.');
end
if (numel(L)~=1)||(numel(L1)~=1)||(numel(L2)~=1)||...
   (numel(M)~=1)||(numel(M1)~=1)||(numel(M2)~=1)
    error('all arguments must have one element.');
end
if (mod(2*L+1,1)~=0)||(mod(2*L1+1,1)~=0)||(mod(2*L2+1,1)~=0)||...
   (mod(2*M+1,1)~=0)||(mod(2*M1+1,1)~=0)||(mod(2*M2+1,1)~=0)
    error('all arguments must be integer or half-integer.');
end
if (L>1e4)||(L1>1e4)||(L2>1e4)
    error('you must be joking.');
end
end

% The miracle of the appropriateness of the language of mathematics for the
% formulation of the laws of physics is a wonderful gift which we neither
% understand nor deserve. We should be grateful for it and hope that it will
% remain valid in future research and that it will extend, for better or for
% worse, to our pleasure, even though perhaps also to our bafflement, to
% wide branches of learning.
%
% Eugene Wigner, http://dx.doi.org/10.1002/cpa.3160130102

