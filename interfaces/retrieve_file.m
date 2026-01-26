% Retrieves a file from an HTTPS link and stores it in a user-
% specified directory. Syntax:
%
%     file_path=retrieve_file(file_url,file_name,dest_dir)
%
% Parameters:
%
%     file_url  - HTTPS URL pointing to the file to retrieve
%
%     file_name - name of the file to be stored on disk
%
%     dest_dir  - destination directory for the downloaded file
%
% Outputs:
%
%     file_path - full path to the retrieved file
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=retrieve_file.m>

function retrieve_file(file_url,file_name,dest_dir)

% Check consistency
grumble(file_url,file_name,dest_dir);

% Ensure destination directory exists
if ~exist(dest_dir,'dir')

    % Create destination directory
    mkdir(dest_dir);

end

% Build destination path
file_path=fullfile(dest_dir,file_name);

% Retrieve the file
websave(file_path,file_url);

end

% Consistency enforcement
function grumble(file_url,file_name,dest_dir)

% Check URL string type
if (~ischar(file_url))&&(~(isstring(file_url)&&isscalar(file_url)))
    error('file_url must be a character string.');
end

% Check file name string type
if (~ischar(file_name))&&(~(isstring(file_name)&&isscalar(file_name)))
    error('file_name must be a character string.');
end

% Check destination string type
if (~ischar(dest_dir))&&(~(isstring(dest_dir)&&isscalar(dest_dir)))
    error('dest_dir must be a character string.');
end

% Check URL scheme
file_url=char(file_url);
if ~strncmp(file_url,'https://',8)
    error('file_url must begin with https://');
end

% Check file name content
file_name=char(file_name);
if isempty(file_name)
    error('file_name must be a non-empty string.');
end

% Check destination content
dest_dir=char(dest_dir);
if isempty(dest_dir)
    error('dest_dir must be a non-empty string.');
end

end

% The "Meir Lahav Criterion" for faculty recruitment at the 
% Weizmann Institute of Science: the applicant must be able
% to change direction without losing momentum.

