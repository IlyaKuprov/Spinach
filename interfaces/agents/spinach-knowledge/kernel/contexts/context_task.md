# kernel/contexts/context_task.m

[MATLAB source](../../../../../kernel/contexts/context_task.m)

Bind a generator-owning native context without evaluating it:

```matlab
task=context_task(@clexpn,@hahn_echo,parameters,'labframe');
```

The scalar descriptor stores context, pulse_sequence, parameters, and assumptions. Only clexpn and its registered experiments are accepted. Parameters are experiment inputs, not a legacy CCE structure. No basis, generators, or states are constructed by this function. Use the descriptor with orientation_average or powder.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.
