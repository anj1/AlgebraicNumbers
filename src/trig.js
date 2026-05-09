import { AlgebraicNumber, root } from './algebraic.js';
import { C } from './complex.js';
import { addPoly, divPolyExact, mulPoly, normalizeAscending, samePolynomialUpToScale, subPoly } from './poly.js';
import { divisors, invTotient } from './invTotient.js';

function gcd(a, b) {
  a = Math.abs(Math.trunc(a));
  b = Math.abs(Math.trunc(b));
  while (b) [a, b] = [b, a % b];
  return a || 1;
}

function normalizeRational(num, den = 1) {
  if (!Number.isInteger(num) || !Number.isInteger(den) || den === 0) throw new RangeError('expected integer numerator and denominator');
  if (den < 0) { num = -num; den = -den; }
  const g = gcd(num, den);
  return [num / g, den / g];
}

export function cyclotomicPolynomial(n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError('n must be a positive integer');
  if (n === 1) return [-1n, 1n];
  let p = Array(n + 1).fill(0n);
  p[0] = -1n;
  p[n] = 1n;
  for (const d of divisors(n).filter(d => d < n)) {
    const { quotient, remainder } = divPolyExact(p, cyclotomicPolynomial(d), 1e-7);
    if (remainder.some(c => c !== 0n)) throw new Error(`failed cyclotomic division for n=${n}, d=${d}`);
    p = quotient;
  }
  return normalizeAscending(p);
}

export function expAlg(num, den = 1) {
  [num, den] = normalizeRational(num, den);

  // exp(pi*i*num/den) = exp(2*pi*i*(num/(2*den))).
  const [, order] = normalizeRational(num, 2 * den);
  const coeff = cyclotomicPolynomial(order);
  const q = num / den;
  return new AlgebraicNumber(coeff, C(Math.cos(Math.PI * q), Math.sin(Math.PI * q)));
}

export function cosAlg(num, den = 1) {
  return expAlg(num, den).realPart();
}

export function sinAlg(num, den = 1) {
  return expAlg(num, den).imagPart();
}

export function isCyclotomicPolynomial(poly) {
  const deg = poly.length - 1;
  for (const n of invTotient(deg)) {
    if (samePolynomialUpToScale(poly, cyclotomicPolynomial(n), 1e-7)) return { ok: true, n };
  }
  return { ok: false, n: 0 };
}

export function logAlg(a) {
  a = AlgebraicNumber.from(a);
  const { ok, n } = isCyclotomicPolynomial(a.coeff);
  if (!ok) return null;
  const thetaOverPi = Math.atan2(a.apprx.im, a.apprx.re) / Math.PI;
  return normalizeRational(Math.round(n * thetaOverPi), n);
}

export function acosAlg(x) {
  x = AlgebraicNumber.from(x);
  const one = AlgebraicNumber.one();
  const y = root(one.sub(x.mul(x)), 2);
  const z = x.add(AlgebraicNumber.i().mul(y));
  return logAlg(z);
}

export function asinAlg(x) {
  x = AlgebraicNumber.from(x);
  const one = AlgebraicNumber.one();
  const y = root(one.sub(x.mul(x)), 2);
  const z = y.add(AlgebraicNumber.i().mul(x));
  return logAlg(z);
}
