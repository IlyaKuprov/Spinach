# kernel/utilities/md5_hash.m

- Signature: `hashstr=md5_hash(A)`

## Purpose

MD5 hash of any Matlab object as a hex string. Identical sparse and full matrices return different hashes. Syntax: hashstr=md5_hash(A)

## Physical / mathematical content

- Produces a hexadecimal MD5 digest of serialized MATLAB object bytes. Objects with different serializations, including sparse and full matrices, can have different hashes.

## Numerical / algorithmic content

- Computes MD5 over the byte stream returned by `serializeToBytes` and renders each digest byte as two hexadecimal digits.

## Parameters / inputs

- A -Matlab object of any type

## Outputs

- hashstr -hexadecimal character string

## Implementation structure

- Serializes `A`, computes its MD5 digest, and converts the digest bytes to a hexadecimal character string.
