# kernel/includes/end_disallow_gpu.m

- Signature: `(script file)`

## Purpose

Reinstates GPU arithmetic setting to the previous state after the start_disallow_gpu command had been issued.

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.

## Numerical / algorithmic content

## Implementation structure

- Reinstates GPU arithmetic setting to the previous state after
- the start_disallow_gpu command had been issued.
- Check that a disallow command had been called
- Return GPU policy to its previous state
- Юлий Ким, "Истерическая
- перестроечная", 1988
- Ну ребята, всё ребята,
- Нету хода нам назад,
- Оборвалися канаты,
- Тормоза не тормозят.
- Вышла фига из кармана,
- Тут же рухнули мосты,
