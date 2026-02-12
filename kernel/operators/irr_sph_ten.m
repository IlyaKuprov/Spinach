% Single-spin irreducible spherical tensor operators T(k,m)  
% obeying the following commutation relation:
%
%                      [Lz,T_km]=m*T_km
% 
% Syntax:
%
%                    T=irr_sph_ten(mult,k)
%
% Parameters:
%
%     mult - multiplicity of the spin in question
%
%        k - irreducible spherical tensor rank (optional)
%
% Outputs:
%
%        T - a two-argument call returns a cell array of 
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
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=irr_sph_ten.m>

function T=irr_sph_ten(mult,k)

% Adapt to the input style
switch nargin
    
    case 1
        
        % Generate tensors of all ranks and put them into a cell array
        if isnumeric(mult)&&isscalar(mult)&&(mult>0)&&(mod(mult,1)==0)
            
            % Preallocate the answer
            T=cell(mult^2,1);
            
            % Fill in recursively
            for n=0:(mult-1)
                T((n^2+1):((n+1)^2))=irr_sph_ten(mult,n);
            end
            
        else
            
            % Catch incorrect calls
            error('mult must be a non-negative integer.');
        
        end
        
    case 2
        
        % Generate spherical tensor operators of the specified rank
        if isnumeric(mult)&&isscalar(mult)&&(mod(mult,1)==0)&&...
           isnumeric(k)&&isscalar(k)&&(k>0)&&(k<mult)&&(mod(k,1)==0)
            
            % Get Pauli matrices
            L=pauli(mult);
            
            % Preallocate the cell array
            T=cell(2*k+1,1);
            
            % Get the top state
            T{1}=((-1)^k)*(2^(-k/2))*L.p^k;
            
            % Apply sequential lowering using Racah's commutation rule
            for n=2:(2*k+1)
                q=k-n+2; T{n}=(1/sqrt((k+q)*(k-q+1)))*(L.m*T{n-1}-T{n-1}*L.m);
            end
                
        elseif isnumeric(mult)&&isscalar(mult)&&(mod(mult,1)==0)&&...
               isnumeric(k)&&isscalar(k)&&(k==0)
            
            % For zero rank, return a unit matrix
            T={speye(mult)};
            
        else
            
            % Catch incorrect calls
            error('mult must be a non-negative integer and k must be an integer from [0, mult-1] interval.');
            
        end
        
end

end

% I swear by my life and my love of it that I will never 
% live for the sake of another man, nor ask another man
% to live for mine.
%
% Ayn Rand, "Atlas Shrugged"

