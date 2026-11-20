# ddf_regularised: softened midpoint field

Use this isolated convolution only when comparing with a solver that samples a softened kernel at source and target cell centres. The regularisation radius is a physical length; changing it changes the model. Unlike the cell-integrated cavity field, the diagonal is either explicitly omitted or retained using zero squared direction cosine. Both choices match their BlochDDFGalerkin counterparts under the angular-frequency conversion in the function header. Applying the symbol with ddf_cartesian does not include diffusion or relaxation.
