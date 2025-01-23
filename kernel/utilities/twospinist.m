% Two-spin irreducible spherical tensor operators. Syntax:
%
%   T=twospinist(spin_system,spin_a,spin_b,indices,type)
%
% Parameters:
%
%     spin_a   - number of the first spin
%
%     spin_b   - number of the second spin
%
%     indices  - two-element vector [L,M] containing
%                rank and projection index; L=1,2 are
%                available
%
% In Liouville space, type can be set to:
%
%    'left' - produces left side product superoperator
%
%    'right' - produces right side product superoperator
%
%    'comm' - produces commutation superoperator (default)
%
%    'acomm' - produces anticommutation superoperator
%
% In Hilbert space calculations, the type parameter is ig-
% nored, and the operator itself is always returned.
%
% Outputs:
%
%     T  - irreducible spherical tensor operator
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=twospinist.m>

function T=twospinist(spin_system,spin_a,spin_b,indices,type)

% Check consistency
grumble(spin_system,spin_a,spin_b,indices,type);

% Indexing switches
switch indices(1)
    
    case 1
        
        % Three first rank operators
        switch indices(2)
            
            case +1
                
                T=-(1/2)*(operator(spin_system,{'L+','Lz'},{spin_a spin_b},type,'csc')-...
                          operator(spin_system,{'Lz','L+'},{spin_a spin_b},type,'csc'));
                    
            case  0
                
                T=-(1/8)*(operator(spin_system,{'L+','L-'},{spin_a spin_b},type,'csc')-...
                          operator(spin_system,{'L-','L+'},{spin_a spin_b},type,'csc'));
                      
            case -1
                
                T=-(1/2)*(operator(spin_system,{'L-','Lz'},{spin_a spin_b},type,'csc')-...
                          operator(spin_system,{'Lz','L-'},{spin_a spin_b},type,'csc'));
                      
            otherwise
                
                % Complain and bomb out
                error('incorrect spherical tensor indices.');
                
        end
        
    case 2
        
        % Five second rank operators
        switch indices(2)
            
            case +2
                
                T=+(1/2)*operator(spin_system,{'L+','L+'},{spin_a spin_b},type,'csc');
                
            case +1
                
                T=-(1/2)*(operator(spin_system,{'Lz','L+'},{spin_a spin_b},type,'csc')+...
                          operator(spin_system,{'L+','Lz'},{spin_a spin_b},type,'csc'));
                      
            case  0
                
                T=+sqrt(2/3)*(operator(spin_system,{'Lz','Lz'},{spin_a spin_b},type,'csc')-...
                       (1/4)*(operator(spin_system,{'L+','L-'},{spin_a spin_b},type,'csc')+...
                              operator(spin_system,{'L-','L+'},{spin_a spin_b},type,'csc')));
                      
            case -1
                
                T=+(1/2)*(operator(spin_system,{'Lz','L-'},{spin_a spin_b},type,'csc')+...
                          operator(spin_system,{'L-','Lz'},{spin_a spin_b},type,'csc'));
                
            case -2
                
                T=+(1/2)*operator(spin_system,{'L-','L-'},{spin_a spin_b},type,'csc');
                
            otherwise
                
                % Complain and bomb out
                error('incorrect spherical tensor indices.');
                
        end
        
    otherwise
        
        % Complain and bomb out
        error('incorrect spherical tensor indices.');
        
end
                
end

% Consistency enforcement
function grumble(spin_system,spin_a,spin_b,indices,type)
if (~isnumeric(spin_a))||(~isscalar(spin_a))||(~isreal(spin_a))||...
   (spin_a<1)||(spin_a>spin_system.comp.nspins)||(mod(spin_a,1)~=0)
    error('spin_a must be a positive integer smaller or equal to the number of spins in the system.');
end
if (~isnumeric(spin_b))||(~isscalar(spin_b))||(~isreal(spin_b))||...
   (spin_b<1)||(spin_b>spin_system.comp.nspins)||(mod(spin_b,1)~=0)
    error('spin_b must be a positive integer smaller or equal to the number of spins in the system.');
end
if spin_a==spin_b
    error('spin_a and spin_b cannot be equal.');
end
if ~ischar(type)
    error('type must be a character string.');
end
if (~isnumeric(indices))||(numel(indices)~=2)||(~isreal(indices))
    error('indices must be a two-element real vector.');
end
if ~ismember(indices(1),[1 2])
    error('spherical ranks other than 1 and 2 are not implemented.');
end
if mod(indices(2),1)~=0
    error('projection index must be an integer.');
end
end

% Oh, please... she certainly is not another Marie Curie, and not Ada
% Lovelace of magnetism. The nearest historical parallel for our new-
% est recruit is Caligula's horse.
%
% Overheard in Magdalen College

