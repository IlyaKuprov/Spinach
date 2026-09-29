# interfaces/retrieve_file.m

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/retrieve_file.m) · [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=retrieve_file.m)

## Interface

`file_path=retrieve_file(file_url,file_name,dest_dir)` downloads the URL to `fullfile(dest_dir,file_name)` and returns that destination path. Each input may be a character array or scalar string; after conversion, the URL must begin literally with `https://`, and `file_name` and `dest_dir` must be non-empty. No other filename or URL-content rule is imposed here.

If `dest_dir` is absent, the function calls `mkdir` before forming the path. The transfer is performed by MATLAB `websave(file_path,file_url)`; download errors are not caught by this wrapper. It returns the full destination path and does not define a separate cache, retry, or unit convention.
