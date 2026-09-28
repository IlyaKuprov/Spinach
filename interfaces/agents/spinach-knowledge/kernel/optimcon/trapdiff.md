# kernel/optimcon/trapdiff.m

- Signature: `[DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)`

Computes derivatives of the trapezium-product interval propagator `expm(-i*((HL+HR)/2 + i*dt*(sqrt(3)/12)*[HL,HR])*dt)` with respect to the left and right control coefficients. The quadrature is described by Iserles and Norsett (Corollary 3.3): https://doi.org/10.1098/rsta.1999.0362. Derivatives use Eq. 16 of Goodwin and Kuprov: https://doi.org/10.1063/1.4928978.

`Hd` contains two same-size square drift generators, at the left and right interval edges; `Hc` is a matching control operator or superoperator. `dt` is the positive interval duration in seconds; `cL` and `cR` are real scalar control coefficients at those edges. `spin_system` is passed to `dirdiff`. Outputs `DL` and `DR` are propagator derivatives with respect to `cL` and `cR`, respectively.