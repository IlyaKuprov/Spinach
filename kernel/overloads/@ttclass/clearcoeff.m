% Absorbs the physical coefficient of the tensor train into
% its cores. Syntax: 
% 
%                     tt=clearcoeff(tt)
%
% The value of the tensor train remains the same.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/clearcoeff.m>

function tt=clearcoeff(tt)

% Get the number of cores and trains
[ncores,ntrains]=size(tt.cores);

% Loop over the trains in the buffer
for n=1:ntrains
    
    % Scale the coefficient
    A=tt.coeff(1,n)^(1/ncores);
    
    % Apply it to cores
    for k=1:ncores
        tt.cores{k,n}=A*tt.cores{k,n};
    end
    
    % Erase the coefficient
    tt.coeff(1,n)=1;
    
end

end

% "Moral outrage is a middle-class luxury."
%
% Anonymous philosophy professor

