# kernel/overloads/@ttclass/numel.m

- Signature: `n=numel(tt)`
- Source: [`kernel/overloads/@ttclass/numel.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/numel.m)
- Wiki: [`ttclass/numel.m`](https://spindynamics.org/wiki/index.php?title=ttclass/numel.m)

## Shape action

Checks that `tt` is a `ttclass`, obtains `sizes(tt)`, converts its dimensions to `int64`, and multiplies all entries of that size array using native integer arithmetic. This counts logical matrix elements represented by the tensor train, not the number of stored core entries; it does not inspect or expand the cores and is independent of coefficient values and ranks.

## Result and guards

If the computed count exceeds MATLAB's `flintmax`, the method errors because that count cannot be represented exactly as a double. Otherwise it converts the integer count to a double scalar. Non-`ttclass` inputs raise an error. No conjugation or core transformation is involved.
