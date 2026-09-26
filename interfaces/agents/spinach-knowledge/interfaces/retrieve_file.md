# interfaces/retrieve_file.m

- Signature: `file_path=retrieve_file(file_url,file_name,dest_dir)`

## Purpose

Retrieves a file from an HTTPS link and stores it in a user- specified directory. Syntax: file_path=retrieve_file(file_url,file_name,dest_dir)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- file_url -HTTPS URL pointing to the file to retrieve
- file_name -name of the file to be stored on disk
- dest_dir -destination directory for the downloaded file

## Outputs

- file_path -full path of the downloaded file on disk

## Implementation structure

- Retrieves a file from an HTTPS link and stores it in a user-
- specified directory. Syntax:
- file_path=retrieve_file(file_url,file_name,dest_dir)
- file_url -HTTPS URL pointing to the file to retrieve
- file_name -name of the file to be stored on disk
- dest_dir -destination directory for the downloaded file
- file_path -full path of the downloaded file on disk
- Check consistency
- Ensure destination directory exists
- Create destination directory
- Build destination path
- Retrieve the file
