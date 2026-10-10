# kernel/overloads/@ttclass/trace.m

## Signature

`tttrace=trace(tt)`

## Behaviour

For every core and every pair of left/right bond indices, the method reshapes that core's physical row/column slice to a matrix of the corresponding local mode sizes and applies MATLAB's `trace` to that matrix. The resulting local traces are stored as cores with both physical modes set to one. The core sequence and bond ranks are retained, as are the input train coefficients; the auxiliary representation's tolerance is set to zero for each train.

Finally, `full` materialises the singleton-mode auxiliary tensor train, yielding the trace as a scalar. The implementation performs these local diagonal contractions and the final train summation directly; it does not invoke rank truncation, rounding, or a tolerance-controlled approximation.

## References

- [Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/trace.m)
- [Spin Dynamics Wiki: `ttclass/trace.m`](https://spindynamics.org/wiki/index.php?title=ttclass/trace.m)
