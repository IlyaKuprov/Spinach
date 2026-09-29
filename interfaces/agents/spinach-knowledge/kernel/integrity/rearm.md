# kernel/integrity/rearm.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/integrity/rearm.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rearm.m)

`rearm()` rebuilds the baseline used by [`sniff`](sniff.md) to flag edits to Spinach MATLAB files. It is an integrity utility, not a physical model.

It recursively scans .m files below `kernel`, `interfaces`, `experiments`, and `etc`. For each included file, it appends a fingerprint formed from the file's basename and a hash of its line contents: `md5_hash([filename md5_hash(content)])`. The current exception list is empty. The fingerprints are stored as the `smells` cell array in `kernel/integrity/smells.mat`; this is a list, not a path-to-hash index. The directory path is not part of a fingerprint.

Call it after establishing the source version that should count as the baseline. A later `sniff` reports fingerprints absent from that list. Because fingerprints omit directory paths and the comparison is membership-only, this scheme does not identify moved files by path, nor does it detect files that have been deleted from the scan tree.

## Inputs and outputs

No inputs or return value. It deletes the existing `smells.mat`, saves the newly collected `smells`, and displays `rearm: sniffer rearmed.`. The source contains no input-validation guard.
