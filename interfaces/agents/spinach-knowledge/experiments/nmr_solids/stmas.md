# Shifted-echo STMAS

`stmas` is a `singlerot` callback for one spin-3/2 87Rb nucleus in `zeeman-liouv`. Three finite x pulses excite satellites, convert to CT+, and invert before CT− detection. The conversion starts at the same rotor phase for each indirect increment; the inversion phase follows the actual elapsed time. Direct samples begin immediately after inversion. The complex result has direct time in rows and indirect time in columns; the dwells are one rotor period and 17/9 rotor periods, respectively. The spectral widths are the reciprocals of those dwells.
