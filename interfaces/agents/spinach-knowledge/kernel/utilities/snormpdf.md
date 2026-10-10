# kernel/utilities/snormpdf.m

## Purpose

Evaluates the probability density of Azzalini's skew normal distribution, given a location, a scale, and a skew factor.

## Behaviour

- Syntax: `p=snormpdf(x,mu,sigma,alpha)`.
- The function first calls an internal consistency checker (`grumble`) on all four inputs.
- The density is computed as `p=2*normpdf(x,mu,sigma).*normcdf(alpha*x,alpha*mu,sigma)`, which the source identifies as Equation 2 in the JSTOR reference.
- `p` has the same shape as `x`.

## Inputs and outputs

Inputs:

- `x` — an array of real numbers; must be a real numeric array.
- `mu` — expectation value of the normal distribution; must be a real scalar.
- `sigma` — standard deviation of the normal distribution; must be a real positive scalar.
- `alpha` — skew factor, a real number; must be a real scalar.

Output:

- `p` — an array of probability densities, same shape as `x`.

Errors are raised with the messages `x must be a real numeric array.`, `mu must be a real scalar.`, `sigma must be a real positive scalar.`, and `alpha must be a real scalar.` when the corresponding checks fail.

## References

- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/snormpdf.m>
- Wiki: <https://spindynamics.org/wiki/index.php?title=snormpdf.m>
- Azzalini skew normal density, Equation 2: <http://www.jstor.org/stable/4615982>
