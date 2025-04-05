% Converts the secular relaxation superoperator into the Redfield
% kite form by dropping all non-longitudinal cross-relaxation pro-
% cesses. Useful when the relaxation superoperator is huge, but
% TROSY-like effects are negligible. Syntax:
%
%                     R=sec2kite(spin_system,R)
%
% Parameters:
%
%    R - relaxation superoperator
%
% Outputs:
%
%    R - relaxation superoperator
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sec2kite.m>

function R=sec2kite(spin_system,R)

% Check consistency
grumble(spin_system,R);

% Get nonzero count
nnz_before=nnz(R);      

% Compile the index of all longitudinal product states in the basis
[~,M]=lin2lm(spin_system.bas.basis); long_states=find(sum(abs(M),2)==0);

% Convert R to XYZ format
[rows,cols,vals]=find(R);

% Zero all rates except self-relaxation and longitudinal cross-relaxation terms
vals=vals.*((ismember(rows,long_states)&ismember(cols,long_states))|(rows==cols));

% Recompose the relaxation superoperator and get nonzero count
R=sparse(rows,cols,vals,length(R),length(R)); nnz_after=nnz(R); 

% Inform the user
report(spin_system,'non-longitudinal cross-relaxation processes dropped,');
report(spin_system,['nnz(R) reduced from ' int2str(nnz_before) ' to ' int2str(nnz_after)]); 

end

% Consistency enforcement
function grumble(spin_system,R)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function requires sphten-liouv formalism.');
end
if (~isnumeric(R))||(size(R,1)~=size(R,2))
    error('R must be a square matrix.');
end
unit=unit_state(spin_system);
if norm(R*unit,2)>1e-10
    error('R appears to be thermalised, cannot proceed.');
end
end

% To Nero, Emperor of Rome, Master of the World, Divine Pontiff. I 
% know that my death will be a disappointment to you, since you wi-
% shed to render me this service yourself. To be born in your reign
% is a miscalculation; but to die in it is a joy. 
%
% I can forgive you for murdering your wife and your mother, for 
% burning our beloved Rome, for befouling our fair country with the
% stench of your crimes. But one thing I cannot forgive - the bore-
% dom of having to listen to your verses, your second-rate songs, 
% your mediocre performances. 
%
% Adhere to your special gifts, Nero - murder and arson, betrayal 
% and terror. Mutilate your subjects if you must; but with my last
% breath I beg you - do not mutilate the arts. Fare well, but com-
% pose no more music. Brutalise the people, but do not bore them,
% as you have bored to death your friend, the late Gaius Petronius.
%
% Petronius, in his dying letter to Nero

