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
% Note: the Schulten-Gordon three-term recurrence is combined with bounded
%       Luscombe-Luban ratios in the nonclassical tails. All intermediate
%       operations use native doubles; coefficients and iterates retain
%       roundoff corrections in a second double. No Java or symbolic toolbox
%       is used. The shortest magnetic-index family is normalised using
%       3j orthogonality. This choice avoids the j=0 recurrence singularity.
%       The existing input ceiling is L,L1,L2<=1e4.
%
% References: https://doi.org/10.1063/1.522426
%             https://doi.org/10.1103/PhysRevE.57.7274
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
            
            % Reject incompatible angular-momentum and projection parities
            if any(mod([a+b+c a+alp b+bet c+gam],1)), return; end

            % Select the shortest 3j family using phase-free cyclic permutations
            spins=[double(a) double(b) double(c)];
            mags=[double(alp) double(bet) -double(gam)];
            lengths=spins([2 3 1])+spins([3 1 2])-...
                    max(abs(spins([2 3 1])-spins([3 1 2])),abs(mags))+1;
            [~,column]=min(lengths);
            order=mod(column-1+(0:2),3)+1;
            spins=spins(order); mags=mags(order);
            mmin=max(-spins(2),-spins(3)-mags(1));
            mmax=min(spins(2),spins(3)-mags(1)); ranks=mmin:mmax;
            count=numel(ranks); target=mags(2)-mmin+1;
            phase=(-1)^(double(a)-double(b)+double(gam)+spins(2)-spins(3)-mags(1));

            % Form the magnetic Schulten-Gordon recurrence coefficients
            edges=mmin:mmax+1;
            factor_a=(spins(2)-edges+1).*(spins(2)+edges);
            factor_b=(spins(3)-edges-mags(1)+1).*(spins(3)+edges+mags(1));
            [product,product_err]=two_product(factor_a,factor_b);
            radial=sqrt(product);
            [squared,squared_err]=two_product(radial,radial);
            radial_err=zeros(size(radial)); nonzero=(radial~=0);
            radial_err(nonzero)=((product(nonzero)-squared(nonzero))+...
                                 (product_err(nonzero)-squared_err(nonzero)))./(2*radial(nonzero));
            upper=radial(2:end); lower=radial(1:end-1);
            centre=spins(2)*(spins(2)+1)+spins(3)*(spins(3)+1)-...
                   spins(1)*(spins(1)+1)-2*ranks.*(ranks+mags(1));

            % Iterate bounded Luscombe-Luban ratios from the left tail
            left_ratio=zeros(2,count); left=1;
            for k=1:count-1
                denom=[centre(k);0];
                if k>1
                    denom=two_add(denom,two_multiply([lower(k);radial_err(k)],left_ratio(:,k-1)));
                end
                if abs(upper(k))>abs(denom(1)), break; end
                left_ratio(:,k)=two_divide([-upper(k);-radial_err(k+1)],denom); left=k+1;
            end

            % Iterate bounded Luscombe-Luban ratios from the right tail
            right_ratio=zeros(2,count); right=count;
            for k=count:-1:2
                denom=[centre(k);0];
                if k<count
                    denom=two_add(denom,two_multiply([upper(k);radial_err(k+1)],right_ratio(:,k+1)));
                end
                if abs(lower(k))>abs(denom(1)), break; end
                right_ratio(:,k)=two_divide([-lower(k);-radial_err(k)],denom); right=k-1;
            end

            % Set the phase without evaluating a potentially underflowing endpoint
            phase=phase*prod(sign(right_ratio(1,right+1:count)));
            values=zeros(2,count);
            if left>=right

                % Join overlapping monotone tails at their common anchor
                left=right; values(1,left)=1;

            else

                % Propagate the left solution through the oscillatory region
                match=floor((left+right)/2);
                values(1,left)=1;
                if left>1
                    values(:,left-1)=left_ratio(:,left-1);
                else
                    values(:,2)=two_divide([-centre(1);0],[upper(1);radial_err(2)]);
                end
                for k=max(left,2):min(match,count-1)
                    values(:,k+1)=recurrence_step(centre(k),values(:,k),...
                                                  [lower(k);radial_err(k)],values(:,k-1),...
                                                  [upper(k);radial_err(k+1)]);
                end

                % Propagate the right solution through the oscillatory region
                backward=zeros(2,count); backward(1,right)=1;
                if right<count
                    backward(:,right+1)=right_ratio(:,right+1);
                else
                    backward(:,count-1)=two_divide([-centre(count);0],[lower(count);radial_err(count)]);
                end
                for k=min(right,count-1):-1:max(match,2)
                    backward(:,k-1)=recurrence_step(centre(k),backward(:,k),...
                                                    [upper(k);radial_err(k+1)],backward(:,k+1),...
                                                    [lower(k);radial_err(k)]);
                end

                % Fit neighbouring values to avoid division at an oscillatory zero
                overlap=max(1,match-1):min(count,match+1);
                scale=sum(values(1,overlap).*backward(1,overlap))/sum(values(1,overlap).^2);
                values(:,left:match)=scale*values(:,left:match);
                values(:,match+1:right)=backward(:,match+1:right);

            end

            % Reconstruct both tails outwards using ratios of magnitude at most one
            for k=left-1:-1:1
                values(:,k)=two_multiply(left_ratio(:,k),values(:,k+1));
            end
            for k=right+1:count
                values(:,k)=two_multiply(right_ratio(:,k),values(:,k-1));
            end

            % Combine weighted 3j normalisation with the CG conversion
            values=sum(values,1);
            cg=phase*values(target)*sqrt((2*double(c)+1)/((2*spins(1)+1)*sum(values.^2)));

        end
            
    end
        
end

end

% Error-free product of native doubles using Dekker splitting
function [product,residual]=two_product(left,right)
product=left.*right;
split_left=(2^27+1)*left; split_right=(2^27+1)*right;
left_hi=split_left-(split_left-left); left_lo=left-left_hi;
right_hi=split_right-(split_right-right); right_lo=right-right_hi;
residual=((left_hi.*right_hi-product)+left_hi.*right_lo+left_lo.*right_hi)+left_lo.*right_lo;
end

% Sum of amplitudes represented by two native doubles
function value=two_add(left,right)
added=left(1)+right(1); recovered=added-left(1);
residual=(left(1)-(added-recovered))+(right(1)-recovered)+left(2)+right(2);
value=[added+residual;residual-((added+residual)-added)];
end

% Product of amplitudes represented by two native doubles
function value=two_multiply(left,right)
[product,residual]=two_product(left(1),right(1));
residual=residual+left(1)*right(2)+left(2)*right(1);
value=[product+residual;residual-((product+residual)-product)];
end

% Quotient of amplitudes represented by two native doubles
function value=two_divide(numer,denom)
quotient=numer(1)/denom(1);
[product,product_err]=two_product(quotient,denom(1));
correction=((numer(1)-product)+numer(2)-product_err-quotient*denom(2))/denom(1);
value=[quotient+correction;correction-((quotient+correction)-quotient)];
end

% Compensated three-term recurrence with two native doubles per amplitude
function value=recurrence_step(centre,current,lower,previous,upper)
value=two_add(two_multiply([centre;0],current),two_multiply(lower,previous));
value=-two_divide(value,upper);
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


