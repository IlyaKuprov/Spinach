# interfaces/retrieve_file.m

- Signature: `file_path=retrieve_file(file_url,file_name,dest_dir)`

## Purpose

Retrieves a file from an HTTPS URL and stores it in a specified directory.

## Parameters / inputs

- `file_url` — HTTPS URL pointing to the file to retrieve. Must be a character array or scalar string and begin with the literal scheme `https://`.
- `file_name` — name to use for the stored file. Must be a non-empty character array or scalar string.
- `dest_dir` — destination directory for the downloaded file. Must be a non-empty character array or scalar string.

## Outputs

- `file_path` — full path of the downloaded file on disk.

## Implementation structure

The function checks the inputs, creates `dest_dir` if it does not exist, builds `file_path` with `fullfile(dest_dir,file_name)`, and downloads the file using `websave(file_path,file_url)`.