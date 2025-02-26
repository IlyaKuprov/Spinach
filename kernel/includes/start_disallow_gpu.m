% Forces GPU arithmetic to be turned off even if the user had 
% requested it in sys.enable setting; restore previous state 
% using end_disallow_gpu command.
%
% ilya.kuprov@weizmann.ac.il

% Check if GPU is currently enabled
user_wanted_gpu=ismember('gpu',spin_system.sys.enable);
report(spin_system,'WARNING: GPU disallowed by programmer request.');

% Disable GPU if it had been enabled
if user_wanted_gpu
    spin_system.sys.enable=setdiff(spin_system.sys.enable,{'gpu'});
end

% Once at MIT, I saw a frat bro accidentally smear a line all over the
% table. Somehow he still managed to snort the whole thing. He looked
% me square in the eyes and said "same high bro, Stokes theorem".
%
% Internet folklore

