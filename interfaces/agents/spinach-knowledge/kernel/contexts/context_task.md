# kernel/contexts/context_task.m

Bind a generator-owning native context without evaluating it:

```matlab
task=context_task(@clexpn,@hahn_echo,parameters,'labframe');
```

The scalar descriptor stores context, pulse_sequence, parameters and assumptions. Only clexpn and its registered experiments are accepted. Parameters are experiment inputs, not a legacy CCE structure. No basis, generators or states are constructed by this function. Use the descriptor with orientation_average or powder.
