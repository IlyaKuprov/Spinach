# kernel/grids/gaussleg.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/gaussleg.m`
- Signature: `[x,w]=gaussleg(a,b,n)`
- Total lines: 75

## Purpose

Computes Gauss-Legendre points and weights in [a,b] interval with accuracy order n. Syntax: [x,w]=gaussleg(a,b,n)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(a,b,n)`.
- Lines 30-31: Initial guess for the nodes in [-1 1]; implemented by `x=cos((2*(0:n)'+1)*pi/(2*n+2))+(0.27/(n+1))*sin(pi*linspace(-1,1,n+1)'*n/(n+2))`.
- Lines 33-34: Newton-Raphson refinement; implemented by `V=zeros(n+1,n+2); prev_x=0`.
- Lines 44-45: Compute the weights; implemented by `w=(b-a)./((1-x.^2).*dV.^2)*((n+2)/(n+1))^2`.
- Lines 47-48: Map from [-1,1] to [a,b]; implemented by `x=(a*(1-x)+b*(1+x))/2`.
- Lines 50-51: Sort arguments in ascending order; implemented by `[x,idx]=sort(x,'ascend'); w=w(idx)`.

### Control flow inferred from the code

- Line 35: `while` loop over `max(abs(x-prev_x))>eps`.
- Line 37: `for` loop over `k=2:(n+1)`.

### Key state/data transformations

- Lines 31: computes `x` using `x=cos((2*(0:n)'+1)*pi/(2*n+2))+(0.27/(n+1))*sin(pi*linspace(-1,1,n+1)'*n/(n+2))`.
- Lines 34: computes `V` using `V=zeros(n+1,n+2); prev_x=0`.
- Lines 36: computes `V(:,1)` using `V(:,1)=1; V(:,2)=x`.
- Lines 38: computes `V(:,k+1)` using `V(:,k+1)=((2*k-1)*x.*V(:,k)-(k-1)*V(:,k-1))/k`.
- Lines 40: computes `dV` using `dV=(n+2)*(V(:,n+1)-x.*V(:,n+2))./(1-x.^2)`.
- Lines 41: computes `prev_x` using `prev_x=x; x=x-V(:,n+2)./dV`.
- Lines 45: computes `w` using `w=(b-a)./((1-x.^2).*dV.^2)*((n+2)/(n+1))^2`.
- Lines 51: computes `[x,idx]` using `[x,idx]=sort(x,'ascend'); w=w(idx)`.

### Local helper functions

- Line 56: `grumble()` — `function grumble(a,b,n)`.
  - Representative operation: `if (~isnumeric(a))||(~isreal(a))||(numel(a)~=1)||(~isfinite(a))`.
  - Representative operation: `error('a must be a finite real number.')`.

## Parameters / inputs

- a -left edge of the interval
- b -right edge of the interval
- n -accuracy order, the number of points in the
- resulting grid will be n+1.

## Outputs

- x -Gauss-Legendre points
- w -Gauss-Legendre weights

## Implementation structure

- Computes Gauss-Legendre points and weights in [a,b] interval
- with accuracy order n. Syntax:
- [x,w]=gaussleg(a,b,n)
- a -left edge of the interval
- b -right edge of the interval
- n -accuracy order, the number of points in the
- resulting grid will be n+1.
- x -Gauss-Legendre points
- w -Gauss-Legendre weights
- Check consistency
- Initial guess for the nodes in [-1 1]
- Newton-Raphson refinement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
