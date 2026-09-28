# kernel/utilities/wigner_6j.m

- Signature: `w=wigner_6j(j1,j2,j3,j4,j5,j6)`

Computes the Wigner 6j-symbol with indices arranged as

```text
/ j1 j2 j3 \
\ j4 j5 j6 /
```

`j1`–`j6` are angular-momentum indices (integers or half-integers); `w` is the resulting 6j-symbol. Physically inadmissible indices give zero.

## Computation

The function initializes `w=0` and sums over `m1=-j1:j1`, …, `m6=-j6:j6`. Only terms satisfying all four screening conditions contribute:

```text
m1+m2-m3=0,  -m1+m5+m6=0,
m4-m5+m3=0,  -m4-m2-m6=0.
```

For each such term, it adds

```text
(-1)^p * wigner_3j(j1,m1,j2,m2,j3,-m3)
       * wigner_3j(j1,-m1,j5,m5,j6,m6)
       * wigner_3j(j4,m4,j5,-m5,j3,m3)
       * wigner_3j(j4,-m4,j2,-m2,j6,-m6)
```

where `p=(j1-m1)+(j2-m2)+(j3-m3)+(j4-m4)+(j5-m5)+(j6-m6)`.

Before summation, `grumble` checks that every argument is numeric, real, scalar, and integer or half-integer; otherwise it raises an error.

Source: <https://spindynamics.org/wiki/index.php?title=wigner_6j.m>