# tmSPICE N-to-CA transfer replay

This reduced example couples the four-spin N–CA–CO–CB model to two sampled spatial RF trajectories. The receive magnitude is already included in the source spatial quadrature weights. The orientation, starting rotor phase, and spatial sample jointly identify each drift/map pair. The maps are interpolated periodically at propagation midpoints; neither the transmitter phase nor the receive magnitude is inferred from a fitted pulse.

The replay checks the unity-map limit and a directional derivative before a short optimisation. It reports signed transfer and the receive-weighted mean integral of squared angular RF amplitude. It is not a reproduction of the published full-grid optimisation.
