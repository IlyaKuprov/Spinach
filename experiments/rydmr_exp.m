% Singlet-singlet RYDMR experiment with exponential recombination
% function (http://dx.doi.org/10.1080/00268979809483134). Syntax:
%
%             A=rydmr_exp(spin_system,parameters,H,R,K)
%
% where H is the Hamiltonian commutation superoperator in zero ex-
% ternal field, R is the relaxation superoperator and K is the che-
% mical kinetics superoperator. The following parameters are requi-
% red:
%
%    parameters.fields  -  row vector of field values, Tesla; the
%                          primary magnet field should be set to
%                          sys.magnet=1 for normalisation purposes
%
%    parameters.rates  -   row vector of singlet recombination
%                          rate constants, Hz
%
%    parameters.electrons - numbers identifying the two electrons
%                           in the isotope list, e.g. [1 2]
%
%    parameters.needs  -  must contain 'zeeman_op', this is an 
%                         instruction to the kernel to provide a
%                         separate Zeeman operator for field sweep 
%                         purposes
%
% Output:
%
%       A - a matrix of singlet yields with dimensions
%           matching the sizes of parameters.rates and
%           parameters.fields
%
% Note: exponential recombination kinetics is built into this func-
%       tion, do not combine with inter.chem.rp_rates parameter.
% 
% i.kuprov@soton.ac.uk
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rydmr_exp.m>

function answer=rydmr_exp(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Isolate the Zeeman operator
Z=parameters.hzeeman; H=H-Z;

% Get the two-electron singlet state
S=singlet(spin_system,parameters.electrons(1),...
                      parameters.electrons(2));

% Localize input arrays
rates=parameters.rates; fields=parameters.fields;

% Preallocate result arrays
N=numel(rates); M=numel(fields);
answer=zeros([N*M 1],'like',1i);

% Decide the formalism
switch spin_system.bas.formalism
    
    case {'sphten-liouv','zeeman-liouv'}

        % Compose Liouvillian
        L=H+1i*R+1i*K;

        % Normalise the singlet
        S=S/norm(S,2);

        % Merged parfor loop
        parfor nm=1:(N*M)
            
            % Extract indices
            [n,m]=ind2sub([N M],nm);
            
            % Assemble Liouvillian
            L_current=L+fields(m)*Z-1i*rates(n)*speye(size(L)); %#ok<PFBNS>
            
            % Compute RYDMR
            answer(nm)=rates(n)*evolution(spin_system,L_current,S,S,[],[],'total');
            
        end
        
    case {'zeeman-hilb'}
        
        % Normalize the singlet
        S=S/norm(S,'fro');
        
        % Get the unit operator
        Id=unit_state(spin_system);
        
        % Merged parfor loop
        parfor nm=1:(N*M)
            
            % Extract indices
            [n,m]=ind2sub([N M],nm);
            
            % Assemble and tidy up Hamiltonian
            H_curr=H+fields(m)*Z; %#ok<PFBNS>
            H_curr=(H_curr+H_curr')/2;
            
            % Compute integration endpoint
            t_end=10/rates(n); %#ok<PFBNS>
            
            % Compute RYDMR (rates inside for roundoff reasons)
            answer(nm)=hdot(S,expmint(spin_system,H_curr,S*rates(n),H_curr+1i*rates(n)*Id,t_end));
            
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

% Return the answer
answer=reshape(real(answer),N,M);
    
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if spin_system.inter.magnet~=1
    error('unit magnet specification (sys.magnet=1) is required for field sweep experiments.');
end
if ~ismember('zeeman_op',parameters.needs)
    error('this function requires a separate Zeeman operator, add ''zeeman_op'' to parameters.needs'); 
end
if ~isfield(parameters,'fields')
    error('magnetic field values must be supplied in parameters.fields variable.');
end
if (~isnumeric(parameters.fields))||(~isreal(parameters.fields))
    error('parameters.fields must be a vector of real numbers.');
end
if ~isfield(parameters,'rates')
    error('singlet recombination rates must be supplied in parameters.rates variable.');
end
if (~isnumeric(parameters.rates))||(~isreal(parameters.rates))||any(parameters.rates<0)
    error('parameters.rates must be a vector of positive real numbers.');
end
if ~isfield(parameters,'electrons')
    error('electrons must be identified in parameters.electrons variable.');
end
if (~isnumeric(parameters.electrons))||(~isreal(parameters.electrons))||...
   any(parameters.electrons<1)||(numel(parameters.electrons)~=2)||...
   any(mod(parameters.electrons,1)~=0)
    error('parameters.electrons must contain two positive integers.');
end
if any(parameters.electrons>spin_system.comp.nspins)
    error('an electron number supplied is greater than the number of spins in the system.');
end
if any(cellfun(@(x)x(1),spin_system.comp.isotopes(parameters.electrons))~='E')
    error('a spin appearing on parameters.electrons is not an electron.');
end
end

% A common mistake that people make when trying to design something
% completely foolproof is to underestimate the ingenuity of complete
% fools.
%
% Douglas Adams

