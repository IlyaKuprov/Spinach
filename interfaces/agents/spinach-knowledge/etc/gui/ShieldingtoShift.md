# etc/gui/ShieldingtoShift.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/ShieldingtoShift.cpp)

The shielding-reference dialog converts a chemical shielding tensor to a chemical shift using the reference for its element. In matrix form the conversion is δ = σ_ref I − σ: subtract each diagonal from the reference and negate each off-diagonal. Do not subtract the reference from every component, or apply the conversion a second time to an interaction already marked as a shift.

For example a shielding diagonal of (180, 185, 190) ppm and a 190 ppm reference gives (10, 5, 0) ppm. Record which reference was used when comparing the resulting shifts to a calculation or measurement; coordinates alone do not determine it.
