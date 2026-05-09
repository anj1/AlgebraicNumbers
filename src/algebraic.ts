import { Complex, C, I, cdist, ComplexInput } from './complex.js';
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
  trimAscending,
  PolyCoeffs
} from './poly.js';
import { factorPolynomial, rationalApprox } from './factor.js';

const EPS = 1e-8;

export interface AlgebraicOpts {
  eps?: number;
  prec?: number;
  simplify?: boolean;
}

export class AlgebraicNumber {
  coeff: PolyCoeffs;
  apprx: Complex;
  prec: number;

  constructor(coeff: PolyCoeffs, apprx: ComplexInput, opts: AlgebraicOpts = {}) {
    this.coeff = normalizeAscending(coeff, opts.eps ?? EPS);
    this.apprx = Complex.from(apprx);
    this.prec = opts.prec ?? calcPrecision(this.coeff);

    if (opts.simplify !== false) {
      const s = simplify(this);
      this.coeff = s.coeff;
      this.apprx = s.apprx;
      this.prec = s.prec;
    }
  }

  static raw(coeff: PolyCoeffs, apprx: ComplexInput, prec?: number): AlgebraicNumber {
    return new AlgebraicNumber(coeff, apprx, { simplify: false, prec });
  }

  static from(x: ComplexInput | AlgebraicNumber): AlgebraicNumber {
    if (x instanceof AlgebraicNumber) return x;
    if (x instanceof Complex || (typeof x === 'object' && x !== null && 're' in x && 'im' in x)) {
      const z = Complex.from(x);
      if (Math.abs(z.im) < EPS) return AlgebraicNumber.from(z.re);
      return AlgebraicNumber.from(z.re).add(AlgebraicNumber.from(z.im).mul(AlgebraicNumber.i()));
    }
    if (typeof x === 'number') {
       const r = rationalApprox(x, 1e8, 1e-10) || { num: Math.round(x * 1e8), den: 1e8 };
       return new AlgebraicNumber([-BigInt(r.num), BigInt(r.den)], C(x, 0));
    }
    if (typeof x === 'bigint') return new AlgebraicNumber([-x, 1n], C(Number(x), 0));
    throw new TypeError(`Cannot convert ${x} to AlgebraicNumber`);
  }

  static rational(num: number | bigint, den: number | bigint = 1): AlgebraicNumber {
    return new AlgebraicNumber([-BigInt(num), BigInt(den)], C(Number(num) / Number(den), 0));
  }

  static zero(): AlgebraicNumber { return AlgebraicNumber.raw([0n, 1n], C(0, 0), 1); }
  static one(): AlgebraicNumber { return AlgebraicNumber.raw([-1n, 1n], C(1, 0), 1); }
  static i(): AlgebraicNumber { return root(AlgebraicNumber.from(-1), 2); }

  degree(): number { return degree(this.coeff); }
  valueOf(): number { return this.apprx.im === 0 ? this.apprx.re : NaN; }
  toComplex(): Complex { return this.apprx.clone(); }
  toString(): string { return this.apprx.toString(); }

  equals(other: ComplexInput | AlgebraicNumber, eps = EPS): boolean {
    const otherAn = AlgebraicNumber.from(other);
    if (!samePolynomialUpToScale(this.coeff, otherAn.coeff, eps)) return false;
    const sep = Math.min(calcPrecision(this.coeff), calcPrecision(otherAn.coeff));
    return cdist(this.apprx, otherAn.apprx) <= Math.max(eps, sep - eps);
  }

  neg(): AlgebraicNumber {
    const coeff = (this.coeff as any[]).map((c, i) => (i % 2 === 0 ? -c : c)) as PolyCoeffs;
    return new AlgebraicNumber(coeff, this.apprx.neg(), { prec: this.prec });
  }

  inv(): AlgebraicNumber {
    if (this.apprx.abs() < EPS) throw new RangeError('division by zero');
    return new AlgebraicNumber(this.coeff.slice().reverse(), this.apprx.inv());
  }

  add(other: ComplexInput | AlgebraicNumber): AlgebraicNumber {
    const otherAn = AlgebraicNumber.from(other);
    return new AlgebraicNumber(composedSum(this.coeff, otherAn.coeff), this.apprx.add(otherAn.apprx));
  }

  sub(other: ComplexInput | AlgebraicNumber): AlgebraicNumber { return this.add(AlgebraicNumber.from(other).neg()); }

  mul(other: ComplexInput | AlgebraicNumber): AlgebraicNumber {
    const otherAn = AlgebraicNumber.from(other);
    if (this.apprx.abs() < EPS || otherAn.apprx.abs() < EPS) return AlgebraicNumber.zero();
    return new AlgebraicNumber(composedProduct(this.coeff, otherAn.coeff), this.apprx.mul(otherAn.apprx));
  }

  div(other: ComplexInput | AlgebraicNumber): AlgebraicNumber { return this.mul(AlgebraicNumber.from(other).inv()); }

  conj(): AlgebraicNumber {
    return new AlgebraicNumber(this.coeff, this.apprx.conj(), { prec: this.prec });
  }

  abs(): AlgebraicNumber {
    return root(this.mul(this.conj()), 2);
  }

  realPart(): AlgebraicNumber {
    return this.add(this.conj()).mul(AlgebraicNumber.rational(1, 2));
  }

  imagPart(): AlgebraicNumber {
    return this.sub(this.conj()).mul(AlgebraicNumber.raw([1n, 0n, 4n], C(0, -0.5), 0.5));
  }

  root(n: number): AlgebraicNumber { return root(this, n); }
  sqrt(): AlgebraicNumber { return root(this, 2); }
  cbrt(): AlgebraicNumber { return root(this, 3); }

  pow2(): AlgebraicNumber { return pow2(this); }

  residual(): Complex { return evalPolyAscending(this.coeff, this.apprx); }
}

export function calcPrecision(coeff: PolyCoeffs): number {
  coeff = trimAscending(coeff);
  if (degree(coeff) <= 1) return Infinity;
  return 0.5 * minRootSeparation(coeff);
}

export function simplify(an: AlgebraicNumber): AlgebraicNumber {
  const coeff = trimAscending(an.coeff);
  if (degree(coeff) <= 1) return AlgebraicNumber.raw(coeff, an.apprx, an.prec);

  const factors = factorPolynomial(coeff).filter(f => degree(f) >= 1);
  if (factors.length <= 1) return AlgebraicNumber.raw(normalizeAscending(coeff), an.apprx, an.prec);

  let best = null;
  for (const factor of factors) {
    const roots = rootsAscending(factor);
    const dist = roots.reduce((m, r) => Math.min(m, cdist(an.apprx, r)), Infinity);
    if (!best || dist < best.dist) best = { factor: factor as PolyCoeffs, dist };
  }

  const f = normalizeAscending(best!.factor);
  return AlgebraicNumber.raw(f, an.apprx, calcPrecision(f));
}

export function root(x: ComplexInput | AlgebraicNumber, n: number): AlgebraicNumber {
  const an = AlgebraicNumber.from(x);
  if (!Number.isInteger(n) || n === 0) throw new RangeError('n must be a non-zero integer');
  if (n === 1) return an;
  if (n < 0) return root(an.inv(), -n);
  const coeff = interleaveZerosAscending(an.coeff, n - 1);
  return new AlgebraicNumber(coeff, an.apprx.powReal(1 / n));
}

export function sqrt(x: ComplexInput | AlgebraicNumber): AlgebraicNumber { return root(x, 2); }
export function cbrt(x: ComplexInput | AlgebraicNumber): AlgebraicNumber { return root(x, 3); }

export function pow2(x: ComplexInput | AlgebraicNumber): AlgebraicNumber {
  const an = AlgebraicNumber.from(x);
  const cfs = an.coeff;
  let p2: PolyCoeffs;
  const isBig = typeof cfs[0] === 'bigint';
  const zero = isBig ? 0n : 0;
  if (cfs.slice(1).filter((_: any, idx: number) => idx % 2 === 0).every((c: any) => c === zero || Math.abs(Number(c)) < EPS)) {
    // p(x) = q(x^2), so q is obtained by taking the even-power coefficients.
    p2 = cfs.filter((_: any, i: number) => i % 2 === 0) as PolyCoeffs;
  } else {
    // Numeric Float64 equivalent of eliminating x from y = x^2 and p(x) = 0.
    let f = polyFromRoots(rootsAscending(cfs).map(r => r.mul(r)));
    if (isBig) {
       const lc1 = BigInt(cfs[cfs.length - 1] as any);
       const lc = lc1 ** BigInt(degree(cfs));
       const numLc = Number(lc);
       p2 = normalizeAscending(f.map(c => BigInt(Math.round(c * numLc)))) as bigint[];
    } else {
       p2 = f as number[];
    }
  }
  return new AlgebraicNumber(p2, an.apprx.mul(an.apprx));
}

export function algRoots(coeff: PolyCoeffs): AlgebraicNumber[] {
  coeff = normalizeAscending(coeff);
  return rootsAscending(coeff).map(r => new AlgebraicNumber(coeff, r));
}

export function confirmAlgebraicNumber(an: ComplexInput | AlgebraicNumber): Complex {
  const anObj = AlgebraicNumber.from(an);
  return evalPolyAscending(anObj.coeff, anObj.apprx);
}

export const zero = (): AlgebraicNumber => AlgebraicNumber.zero();
export const one = (): AlgebraicNumber => AlgebraicNumber.one();
export const add = (a: ComplexInput | AlgebraicNumber, b: ComplexInput | AlgebraicNumber): AlgebraicNumber => AlgebraicNumber.from(a).add(b);
export const sub = (a: ComplexInput | AlgebraicNumber, b: ComplexInput | AlgebraicNumber): AlgebraicNumber => AlgebraicNumber.from(a).sub(b);
export const mul = (a: ComplexInput | AlgebraicNumber, b: ComplexInput | AlgebraicNumber): AlgebraicNumber => AlgebraicNumber.from(a).mul(b);
export const div = (a: ComplexInput | AlgebraicNumber, b: ComplexInput | AlgebraicNumber): AlgebraicNumber => AlgebraicNumber.from(a).div(b);
export const conj = (a: ComplexInput | AlgebraicNumber): AlgebraicNumber => AlgebraicNumber.from(a).conj();
export const real = (a: ComplexInput | AlgebraicNumber): AlgebraicNumber => AlgebraicNumber.from(a).realPart();
export const imag = (a: ComplexInput | AlgebraicNumber): AlgebraicNumber => AlgebraicNumber.from(a).imagPart();
