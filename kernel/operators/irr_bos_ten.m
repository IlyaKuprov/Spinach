% Single-mode bosonic monomial operators:
%
%                B(k,q) = (Cr^(q))*(An^(k-q))
%
% with q=0:k, obeying the following commutation relations with
% the population number operator N:
%
%                          [specify]
% 
% Syntax:
%
%                   B=irr_sph_ten(nlevels,k)
%
% Parameters:
%
%     nlevels - number of bosonic ladder 
%               population levels
%
%           k - total monomial power (optional)
%
% Outputs:
%
%        B - a two-argument call returns a cell array of 
%            monomials of total power k in the order of
%            increasing q. A single argument call produces
%            monomials of all ranks and puts them into a 
%            cell array in the order of increasing rank,
%            and increasing q within each rank.
%
% Note: operator normalisation in spin dynamics is not a good
% idea. The only way to make the formalism independent of the
% spin quantum number is to impose identical commutation rela-
% tions rather than equal matrix norms.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=irr_bos_ten.m>

function B=irr_bos_ten(nlevels,k)

% Adapt to the input style
switch nargin
    
    case 1
        
        % Generate monomials of all total powers and put them into a cell array
        if isnumeric(nlevels)&&isscalar(nlevels)&&(nlevels>0)&&(mod(nlevels,1)==0)
            
            % Preallocate the answer
            B=cell(nlevels*(nlevels+1)/2,1);
            
            % Fill in recursively
            for n=0:(nlevels-1)
                lower_idx=1+n*(n+1)/2; upper_idx=(n+1)*(n+2)/2;
                B(lower_idx:upper_idx)=irr_bos_ten(nlevels,n);
            end
            
        else
            
            % Catch incorrect calls
            error('nlevels must be a non-negative integer.');
        
        end
        
    case 2
       
        % Generate monomials of the specified total power
        if isnumeric(nlevels)&&isscalar(nlevels)&&(mod(nlevels,1)==0)&&...
           isnumeric(k)&&isscalar(k)&&(k>0)&&(k<nlevels)&&(mod(k,1)==0)
            
            % Get Weyl matrices
            A=weyl(nlevels);
            
            % Preallocate cell array
            B=cell(k+1,1);
            
            % Make the monomial
            for q=0:k
                B{q+1}=(A.c^(q))*(A.a^(k-q));
            end
                
        elseif isnumeric(nlevels)&&isscalar(nlevels)&&(mod(nlevels,1)==0)&&...
               isnumeric(k)&&isscalar(k)&&(k==0)
            
            % For zero power, return a unit matrix
            B={speye(nlevels)};
            
        else
            
            % Catch incorrect calls
            error('nlevels must be a non-negative integer and k must be an integer from [0, nlevels-1] interval.');
            
        end
        
end

end

% "Oh, so you really exist. I thought Littlewood 
%  was a pseudonym Hardy used for his less impor-
%  tant articles."
%
% Norbert Wiener, to John Littlewood

