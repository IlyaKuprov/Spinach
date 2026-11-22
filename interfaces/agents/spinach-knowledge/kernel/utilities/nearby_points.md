# kernel/utilities/nearby_points.m

Cutoff-sized spatial bins restrict distance checks to neighbouring occupied cells. The result indexes both point arrays and includes coincident points and the cutoff shell; callers impose their own self-image and strict-cutoff exclusions.
