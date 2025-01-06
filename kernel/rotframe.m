% Rotating frame transformation with respect to specified spins
% to specified order in perturbation theory, using the formalism
% described in https://doi.org/10.1063/1.4928978 Syntax:
%
%        H=rotframe(spin_system,H0,H,isotope,order)
%
% Parameters:
%
%    H0     - carrier Hamiltonian with respect to which the
%             rotating frame transformation is to be done
%
%    H      - laboratory frame Hamiltonian H0+H1 that is to
%             be transformed into the rotating frame 
%
%    isotope - string, such as '1H', specifying the spins
%              with respect to which the transformation is
%              being computed
%
%    order   - perturbation theory order in the rotating
%              frame transformation, this may be inf
%
% Outputs:
%
%    H       - rotating frame Hamiltonian
%
% Notes: the auxiliary matrix method is massively faster than
%        either commutator series or diagonalisation.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rotframe.m>

function Hr=rotframe(spin_system,H0,H,isotope,order)

% Check consistency
grumble(spin_system,H0,H,isotope,order);

% Compute the period
switch spin_system.bas.formalism
    
    case {'zeeman-liouv','sphten-liouv'}
        
        % Liouville space period for H0
        T=-2*pi/(spin(isotope)*spin_system.inter.magnet);
        
    case {'zeeman-hilb'}
        
        % Hilbert space period for H0
        T=-4*pi/(spin(isotope)*spin_system.inter.magnet);
        
end

% Run the interaction representation transformation
Hr=intrep(spin_system,H0,H,T,order);

end

% Consistency enforcement
function grumble(spin_system,H0,H,isotope,order)
if ~ischar(isotope)
    error('isotope must be a character string');
end
if ~isfield(spin_system.inter,'assumptions')
    error('assumption information is missing, call assume() first.');
end
if ismember(spin_system.inter.assumptions,{'nmr'})
    error('all spins are already in the rotating frame.');
end
if ismember(spin_system.inter.assumptions,{'esr','deer','deer-zz'})&&...
   (isotope(1)=='E')
    error('electrons are already in the rotating frame.');
end
[~,mult]=spin(isotope);
if ismember(spin_system.inter.assumptions,{'qnmr'})&&(mult==2)
    error('spin 1/2 nuclei are already in the rotating frame.');
end
if (~ishermitian(H))||(~ishermitian(H0))
    error('both H and C must be Hermitian.');
end
if ((~isreal(order))||(order<1)||(mod(order,1)~=0))&&(~isinf(order))
    error('order must be a positive real integer.');
end
end

% Finally, I acknowledge the financial support of EPSRC in its best,
% that is, the responsive mode. This Nobel Prize would be absolutely
% impossible without this mode.  [...]  However, I can offer no nice
% words for the EU Framework Programmes which, except for the Europe-
% an Research Council, can be praised only by europhobes for discre-
% diting the whole idea of an effectively working Europe.
%
% Andre Geim's Nobel Lecture, 2010

