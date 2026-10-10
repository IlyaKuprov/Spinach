# kernel/overloads/@polyadic/gpuArray.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/gpuArray.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/gpuArray.m)

`gpuArray(p)` accepts a polyadic object only; its local consistency check raises an error for a non-polyadic input. It passes numeric core factors, prefixes, and suffixes to MATLAB's `gpuArray`, recursively uploading nested polyadics while leaving function-handle cores unchanged. The enclosing object stays polyadic. It does not form the Kronecker products, sum terms, or materialise the represented matrix, and does not add dimension checks or broadcasting. Numeric multiplier cores are therefore stored on the device after this conversion rather than uploaded by their action handles.

The transformation is component-wise: this method neither calls `full` nor defines a separate sparse-storage conversion. Device support and storage details therefore follow MATLAB's `gpuArray` handling of each component. The source cautions that GPUs may be poor at permuting array dimensions and recommends checking whether CPU execution is faster; it does not promise a speedup.
