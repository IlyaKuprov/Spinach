% Tests Matlab installation for basic parallel pool sanity
% and stability. This can be a problem on NMRBox where the
% user does not have admin rights and the software deploy-
% ment automation frequently messes things up.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=basic_sanity_check.m>

function basic_sanity_check()

% Kill the parallel pool
delete(gcp('nocreate'));

% Test parallel pool infrastructure
parpool(feature('numcores'));
parfor n=1:1000
    A=randn(1000); A^4; %#ok<VUNUS>
end
spmd
    A=randn(10000,'codistributed'); A=A^2;
end
A=gather(A); clear('A'); %#ok<NASGU>
store=gcp('nocreate').ValueStore;
store.put({'ident'},{randn(100)});

% Test Java infrastructure
engine=java.security.MessageDigest.getInstance('MD5');
engine.update(getByteStreamFromArray('message'));
engine.digest;

% Spinach tests
existentials();

end

% Hopeless? It's not hopeless.
% Doubtful, but not hopeless.
%
% A-ha

% #NHEAD #NGRUM