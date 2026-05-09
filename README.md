# algebraic-numbers

A JavaScript port of the Julia `AlgebraicNumber` prototype using `BigInt` for exact polynomials and `Float64` for numerical approximations. It represents an algebraic number as:

- `coeff`: ascending polynomial coefficients represented as `BigInt`s, e.g. `[-2n, 0n, 1n]` for `x^2 - 2`
- `apprx`: a selected complex Float64 approximation of the intended root
- `prec`: half the minimum separation between roots of `coeff`

This intentionally does **not** try to reproduce Nemo/BigFloat exactness. It is a numeric library with symbolic factorization used only to reduce reducible polynomials and keep branch identity stable.

## Install

```bash
npm install
```

Dependencies:

- `flo-poly` for real polynomial root isolation / solving
- `algebraic-numbers` for symbolic polynomial factorization

## Run tests

```bash
npm test
```

## Usage

```js
import { AlgebraicNumber, sqrt, algRoots, cosAlg } from 'algebraic-numbers';

const two = AlgebraicNumber.from(2);
const s2 = sqrt(two);

console.log(s2.coeff);       // [-2n, 0n, 1n], i.e. x^2 - 2
console.log(s2.apprx.re);    // 1.4142135623730951
console.log(s2.mul(s2).equals(two)); // true

const roots = algRoots([-6n, 11n, -6n, 1n]); // x^3 - 6x^2 + 11x - 6
console.log(roots.map(r => r.apprx.re)); // approximately [1, 2, 3]

console.log(cosAlg(1, 3).apprx.re); // cos(pi / 3) = 0.5
```

## API notes

The public polynomial convention follows the Julia code: coefficients are ascending by degree.

```js
[-6n, 11n, -6n, 1n] // -6 + 11x - 6x^2 + x^3
```

`flo-poly` expects descending coefficients internally, so this package converts before calling it.

Arithmetic is available as methods:

```js
x.add(y)
x.sub(y)
x.mul(y)
x.div(y)
x.neg()
x.inv()
x.conj()
x.realPart()
x.imagPart()
x.sqrt()
x.cbrt()
```

Exported functional wrappers are also available: `add`, `sub`, `mul`, `div`, `sqrt`, `cbrt`, `root`, `real`, `imag`, and `conj`.

## Limitations

- Float64 arithmetic is not a substitute for exact algebraic-number arithmetic.
- `flo-poly` focuses on real roots; this package uses a Durand-Kerner fallback for complex roots.
- Factorization is best-effort. algebraic-numbers is used first; a rational-root fallback handles common reducible cases such as `x^2 - 4`.
- High-degree composed sums/products can be ill-conditioned and expensive because arithmetic expands pairwise roots numerically.
