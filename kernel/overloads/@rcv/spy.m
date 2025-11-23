% Plot the sparsity pattern of an RCV matrix.
%
% m.keitel@soton.ac.uk

function spy(A)

    spy(sparse(A));

end

% Downloaded a virus for Linux lately and unpacked it. Tried to run it as
% root, didn't work. Googled for 2 hours, found out that instead of
% /usr/local/bin the virus unpacked to /usr/bin for which the user malware
% doesn't have any write permissions, therefore the virus couldn't create a
% process file. Found patched .configure and .make files on some Chinese
% forum, recompiled and rerun it. The virus said it needs the library
% cmalw-lib-2.0. Turns out cmalw-lib-2.0 is shipped with CentOS but not
% with Ubuntu. Googled for hours again and found an instruction to build a
% .deb package from source. The virus finally started, wrote some logs, made
% a core dump and crashed. After 1 hour of going through the logs I disco-
% vered the virus assumed it was running on ext4 and called into its disk
% encryption API. Under btrfs this API is deprecated. The kernel noticed
% and made this partition read-only. Opened the sources, grep'ed the Bit-
% coin wallet and sent $5 out of pity.
%
% Internet folklore

