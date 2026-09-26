# kernel/utilities/snormpdf.m

- Signature: `p=snormpdf(x,mu,sigma,alpha)`

## Purpose

Azzalini's skew normal distribution. Syntax: p=snormpdf(x,mu,sigma,alpha)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- x -an array of real numbers
- mu -expectation value of the normal distribution
- sigma -standard deviation of the normal distribution
- alpha -skew factor, a real number

## Outputs

- p -an array of probability densities,
- same shape as x

## Implementation structure

- Azzalini's skew normal distribution. Syntax:
- p=snormpdf(x,mu,sigma,alpha)
- x -an array of real numbers
- mu -expectation value of the normal distribution
- sigma -standard deviation of the normal distribution
- alpha -skew factor, a real number
- p -an array of probability densities,
- same shape as x
- Check consistency
- Equation 2 in http://www.jstor.org/stable/4615982
- Consistency enforcement
- The smallest minority on earth is the individual. Those who deny
