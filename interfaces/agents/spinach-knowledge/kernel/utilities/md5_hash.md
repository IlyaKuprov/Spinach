# kernel/utilities/md5_hash.m

## Purpose

Returns an MD5 hash of any MATLAB object as a hexadecimal string, per the header comment of [`kernel/utilities/md5_hash.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/md5_hash.m). The header notes that identical sparse and full matrices return different hashes.

## Behaviour

The function serialises the input object into a bytestream with `serializeToBytes`, computes the MD5 digest with `digestMD5`, and formats the digest bytes as a lowercase hexadecimal string via `sprintf('%.2x',...)`.

## Inputs and outputs

- `A` — MATLAB object of any type.
- `hashstr` — hexadecimal character string.

## References

- [Spinach Wiki: md5_hash.m](https://spindynamics.org/wiki/index.php?title=md5_hash.m)
- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/md5_hash.m)
