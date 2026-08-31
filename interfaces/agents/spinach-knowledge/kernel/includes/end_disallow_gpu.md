# kernel/includes/end_disallow_gpu.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/includes/end_disallow_gpu.m`
- Signature: `(script file)`
- Total lines: 67

## Purpose

Reinstates GPU arithmetic setting to the previous state after the start_disallow_gpu command had been issued.

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 1-9: Reinstates GPU arithmetic setting to the previous state after the start_disallow_gpu command had been issued.; implemented by `if ~exist('user_wanted_gpu','var')`.
- Lines 8-9: Check that a disallow command had been called; implemented by `if ~exist('user_wanted_gpu','var')`.
- Lines 13-14: Return GPU policy to its previous state; implemented by `if user_wanted_gpu`.

### Control flow inferred from the code

- Line 9: conditional branch on `~exist('user_wanted_gpu','var')`.
- Line 14: conditional branch on `user_wanted_gpu`.

### Key state/data transformations

- Lines 15: computes `spin_system.sys.enable{end+1}` using `spin_system.sys.enable{end+1}='gpu'`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`.
