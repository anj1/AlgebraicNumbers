export const EPS = 1e-10;

export class Complex {
  constructor(re = 0, im = 0) {
    this.re = Number(re);
    this.im = Number(im);
  }

  static from(x) {
    if (x instanceof Complex) return x;
    if (typeof x === 'number') return new Complex(x, 0);
    if (typeof x === 'bigint') return new Complex(Number(x), 0);
    if (x && typeof x.re === 'number' && typeof x.im === 'number') return new Complex(x.re, x.im);
    throw new TypeError(`Cannot convert ${x} to Complex`);
  }

  clone() { return new Complex(this.re, this.im); }
  add(z) { z = Complex.from(z); return new Complex(this.re + z.re, this.im + z.im); }
  sub(z) { z = Complex.from(z); return new Complex(this.re - z.re, this.im - z.im); }
  neg() { return new Complex(-this.re, -this.im); }

  mul(z) {
    z = Complex.from(z);
    return new Complex(this.re * z.re - this.im * z.im, this.re * z.im + this.im * z.re);
  }

  div(z) {
    z = Complex.from(z);
    const d = z.re * z.re + z.im * z.im;
    if (d === 0) throw new RangeError('Complex division by zero');
    return new Complex((this.re * z.re + this.im * z.im) / d, (this.im * z.re - this.re * z.im) / d);
  }

  inv() {
    const d = this.re * this.re + this.im * this.im;
    if (d === 0) throw new RangeError('Complex inversion of zero');
    return new Complex(this.re / d, -this.im / d);
  }

  conj() { return new Complex(this.re, -this.im); }
  abs2() { return this.re * this.re + this.im * this.im; }
  abs() { return Math.hypot(this.re, this.im); }
  arg() { return Math.atan2(this.im, this.re); }

  powReal(t) {
    const r = this.abs();
    const theta = this.arg();
    const rt = Math.pow(r, t);
    return new Complex(rt * Math.cos(theta * t), rt * Math.sin(theta * t));
  }

  equals(z, eps = EPS) {
    z = Complex.from(z);
    return this.sub(z).abs() <= eps;
  }

  clean(eps = EPS) {
    const re = Math.abs(this.re) < eps ? 0 : this.re;
    const im = Math.abs(this.im) < eps ? 0 : this.im;
    return new Complex(re, im);
  }

  toString(precision = 15) {
    const z = this.clean();
    const re = Number(z.re.toPrecision(precision));
    const im = Number(z.im.toPrecision(precision));
    if (im === 0) return String(re);
    if (re === 0) return `${im}im`;
    return `${re}${im < 0 ? '' : '+'}${im}im`;
  }
}

export const C = (re = 0, im = 0) => new Complex(re, im);
export const I = new Complex(0, 1);

export function cexpi(theta) {
  return new Complex(Math.cos(theta), Math.sin(theta));
}

export function cdist(a, b) {
  return Complex.from(a).sub(b).abs();
}
