# Snapshot POD on a coupled proton pair

The example trains on five proton-offset trajectories and evaluates four interleaved held-out offsets of a damped strongly coupled pair. It compares explicit ranks three and four against ordinary full Spinach evolution and reports complex FID and identically apodised spectral errors, without phase or amplitude fitting.

Training samples, offset coverage, acquisition, and Euclidean basis are explicit. Rank four spans the accessible coherence sector for this compact model; rank-three truncation illustrates that snapshot residual and physical observable error are different diagnostics. This is not evidence of a production speedup or universal error bound. Slow Schur selection, pulse fitting, and offset interpolation are not part of the example.
