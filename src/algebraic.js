import { Complex, C, I, cdist } from './complex.js';
import {
  cleanRealCoeffs,
  composedProduct,
  composedSum,
  degree,
  evalPolyAscending,
  interleaveZerosAscending,
  minRootSeparation,
  normalizeAscending,
  polyFromRoots,
  rootsAscending,
  samePolynomialUpToScale,
  trimAscending
} from './poly.js';
import { factorPolynomial } from './factor.js';

const EPS = 1e-8;

export class AlgebraicNumber {
  constructor(coeff, apprx, opts = {}) {
    this.coeff = cleanRealCoeffs(trimAscending(coeff), opts.eps ?? EPS);
    this.apprx = Complex.from(apprx);
    this.prec = opts.prec ?? calcPrecision(this.coeff);

    if (opts.simplify !== false) {
      const s = simplify(this);
      this.coeff = s.coeff;
      this.apprx = s.apprx;
      this.prec = s.prec;
    }
  }

  static raw(coeff, apprx, prec = undefined) {
    return new AlgebraicNumber(coeff, apprx, { simplify: false, prec });
  }

  static from(x) {
    if (x instanceof AlgebraicNumber) return x;
    if (x instanceof Complex || (x && typeof x.re === 'number' && typeof x.im === 'number')) {
      const z = Complex.from(x);
      if (Math.abs(z.im) < EPS) return AlgebraicNumber.from(z.re);
      return AlgebraicNumber.from(z.re).add(AlgebraicNumber.from(z.im).mul(AlgebraicNumber.i()));
    }
    if (typeof x === 'number') return new AlgebraicNumber([-x, 1], C(x, 0));
    throw new TypeError(`Cannot convert ${x} to AlgebraicNumber`);
  }

  static rational(num, den = 1) {
    if (!Number.isFinite(num) || !Number.isFinite(den) || den === 0) throw new RangeError('invalid rational');
    return new AlgebraicNumber([-num, den], C(num / den, 0));
  }

  static zero() { return AlgebraicNumber.raw([0, 1], C(0, 0), 1); }
  static one() { return AlgebraicNumber.raw([-1, 1], C(1, 0), 1); }
  static i() { return root(AlgebraicNumber.from(-1), 2); }

  degree() { return degree(this.coeff); }
  valueOf() { return this.apprx.im === 0 ? this.apprx.re : NaN; }
  toComplex() { return this.apprx.clone(); }
  toString() { return this.apprx.toString(); }

  equals(other, eps = EPS) {
    other = AlgebraicNumber.from(other);
    if (!samePolynomialUpToScale(this.coeff, other.coeff, eps)) return false;
    const sep = Math.min(calcPrecision(this.coeff), calcPrecision(other.coeff));
    return cdist(this.apprx, other.apprx) <= Math.max(eps, sep - eps);
  }

  neg() {
    const coeff = this.coeff.map((c, i) => (i % 2 === 0 ? -c : c));
    return new AlgebraicNumber(coeff, this.apprx.neg(), { prec: this.prec });
  }

  inv() {
    if (this.apprx.abs() < EPS) throw new RangeError('division by zero');
    return new AlgebraicNumber(this.coeff.slice().reverse(), this.apprx.inv());
  }

  add(other) {
    other = AlgebraicNumber.from(other);
    return new AlgebraicNumber(composedSum(this.coeff, other.coeff), this.apprx.add(other.apprx));
  }

  sub(other) { return this.add(AlgebraicNumber.from(other).neg()); }

  mul(other) {
    other = AlgebraicNumber.from(other);
    if (this.apprx.abs() < EPS || other.apprx.abs() < EPS) return AlgebraicNumber.zero();
    return new AlgebraicNumber(composedProduct(this.coeff, other.coeff), this.apprx.mul(other.apprx));
  }

  div(other) { return this.mul(AlgebraicNumber.from(other).inv()); }

  conj() {
    return new AlgebraicNumber(this.coeff, this.apprx.conj(), { prec: this.prec });
  }

  abs() {
    return root(this.mul(this.conj()), 2);
  }

  realPart() {
    return this.add(this.conj()).mul(AlgebraicNumber.rational(1, 2));
  }

  imagPart() {
    return this.sub(this.conj()).mul(AlgebraicNumber.raw([1, 0, 4], C(0, -0.5), 0.5));
  }

  root(n) { return root(this, n); }
  sqrt() { return root(this, 2); }
  cbrt() { return root(this, 3); }

  pow2() { return pow2(this); }

  residual() { return evalPolyAscending(this.coeff, this.apprx); }
}

export function calcPrecision(coeff) {
  coeff = trimAscending(coeff);
  if (degree(coeff) <= 1) return Infinity;
  return 0.5 * minRootSeparation(coeff);
}

export function simplify(an) {
  const coeff = trimAscending(an.coeff);
  if (degree(coeff) <= 1) return AlgebraicNumber.raw(coeff, an.apprx, an.prec);

  const factors = factorPolynomial(coeff).filter(f => degree(f) >= 1);
  if (factors.length <= 1) return AlgebraicNumber.raw(normalizeAscending(coeff), an.apprx, an.prec);

  let best = null;
  for (const factor of factors) {
    const roots = rootsAscending(factor);
    const dist = roots.reduce((m, r) => Math.min(m, cdist(an.apprx, r)), Infinity);
    if (!best || dist < best.dist) best = { factor, dist };
  }

  const f = normalizeAscending(best.factor);
  return AlgebraicNumber.raw(f, an.apprx, calcPrecision(f));
}

export function root(x, n) {
  const an = AlgebraicNumber.from(x);
  if (!Number.isInteger(n) || n === 0) throw new RangeError('n must be a non-zero integer');
  if (n === 1) return an;
  if (n < 0) return root(an.inv(), -n);
  const coeff = interleaveZerosAscending(an.coeff, n - 1);
  return new AlgebraicNumber(coeff, an.apprx.powReal(1 / n));
}

export function sqrt(x) { return root(x, 2); }
export function cbrt(x) { return root(x, 3); }

export function pow2(x) {
  const an = AlgebraicNumber.from(x);
  const cfs = an.coeff;
  let p2;
  if (cfs.slice(1).filter((_, idx) => idx % 2 === 0).every(c => Math.abs(c) < EPS)) {
    // p(x) = q(x^2), so q is obtained by taking the even-power coefficients.
    p2 = cfs.filter((_, i) => i % 2 === 0);
  } else {
    // Numeric Float64 equivalent of eliminating x from y = x^2 and p(x) = 0.
    p2 = polyFromRoots(rootsAscending(cfs).map(r => r.mul(r)));
  }
  return new AlgebraicNumber(p2, an.apprx.mul(an.apprx));
}

export function algRoots(coeff) {
  coeff = normalizeAscending(coeff);
  return rootsAscending(coeff).map(r => new AlgebraicNumber(coeff, r));
}

export function confirmAlgebraicNumber(an) {
  an = AlgebraicNumber.from(an);
  return evalPolyAscending(an.coeff, an.apprx);
}

export const zero = () => AlgebraicNumber.zero();
export const one = () => AlgebraicNumber.one();
export const add = (a, b) => AlgebraicNumber.from(a).add(b);
export const sub = (a, b) => AlgebraicNumber.from(a).sub(b);
export const mul = (a, b) => AlgebraicNumber.from(a).mul(b);
export const div = (a, b) => AlgebraicNumber.from(a).div(b);
export const conj = a => AlgebraicNumber.from(a).conj();
export const real = a => AlgebraicNumber.from(a).realPart();
export const imag = a => AlgebraicNumber.from(a).imagPart();
