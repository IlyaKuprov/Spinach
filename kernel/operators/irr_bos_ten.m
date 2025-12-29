% Single-mode irreducible bosonic tensor operators B(k,q)  
% obeying the following commutation relation with the po-
% pulation number operator N:
%
%                      [N,B_kq]=q*B_kq
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
%           k - bosonic tensor rank (optional)
%
% Outputs:
%
%        B - a two-argument call returns a cell array of 
%            tensors of rank k in the order of decreasing
%            projection. A single argument call produces
%            tensors of all ranks and puts them into a 
%            cell array in the order of increasing rank,
%            and decreasing projection within each rank.
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
        
        % Generate tensors of all ranks and put them into a cell array
        if isnumeric(nlevels)&&isscalar(nlevels)&&(nlevels>0)&&(mod(nlevels,1)==0)
            
            % Preallocate the answer
            B=cell(nlevels^2,1);
            
            % Fill in recursively
            for n=0:(nlevels-1)
                B((n^2+1):((n+1)^2))=irr_bos_ten(nlevels,n);
            end
            
        else
            
            % Catch incorrect calls
            error('nlevels must be a non-negative integer.');
        
        end
        
    case 2
        
        % Generate bosonic tensor operators of the specified rank
        if isnumeric(nlevels)&&isscalar(nlevels)&&(mod(nlevels,1)==0)&&...
           isnumeric(k)&&isscalar(k)&&(k>0)&&(k<nlevels)&&(mod(k,1)==0)
            
            % Get Weyl matrices
            A=weyl(nlevels);
            
            % Preallocate cell array
            B=cell(2*k+1,1);
            
            % Do not use commutators
            for q=-k:k
                B{q+k+1}=(A.c^(k+q))*(A.a^(k-q));
            end
                
        elseif isnumeric(nlevels)&&isscalar(nlevels)&&(mod(nlevels,1)==0)&&...
               isnumeric(k)&&isscalar(k)&&(k==0)
            
            % For zero rank, return a unit matrix
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

