# etc/gui/FilterFile.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/FilterFile.cpp)

This import dialog changes the retained model, unlike a slider that merely hides interactions. J couplings survive when the absolute scalar (first diagonal element) is strictly greater than the J threshold. Hyperfine and quadrupolar tensors survive when their Frobenius norm, sqrt(trace(MᵀM)), is at least their respective threshold. Equality is therefore handled differently for the scalar and tensor filters. A tensor with zero trace can still pass the norm filter.

The histograms summarize the input distributions. Removing orphan atoms is optional: the check uses retained physical interactions, deliberately ignoring automatically calculated chemical bonds. A nucleus connected by a drawn bond alone can consequently disappear. Compare atom as well as interaction counts after applying a filter.
