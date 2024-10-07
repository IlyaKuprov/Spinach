% Tests NMRBox Matalb installation for basic parallel
% pool sanity and stability.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=nmrbox_sanity_check.m>

function nmrbox_sanity_check()

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