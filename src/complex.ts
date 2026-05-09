export const EPS = 1e-10;

export type ComplexInput = Complex | { re: number, im: number } | number | bigint;

export class Complex {
  re: number;
  im: number;

  constructor(re: number | string = 0, im: number | string = 0) {
    this.re = Number(re);
    this.im = Number(im);
  }

  static from(x: ComplexInput): Complex {
    if (x instanceof Complex) return x;
    if (typeof x === 'number') return new Complex(x, 0);
    if (typeof x === 'bigint') return new Complex(Number(x), 0);
    if (typeof x === 'object' && x !== null && 're' in x && 'im' in x && typeof x.re === 'number' && typeof x.im === 'number') return new Complex(x.re, x.im);
    throw new TypeError(`Cannot convert ${x} to Complex`);
  }

  clone(): Complex { return new Complex(this.re, this.im); }
  add(z: ComplexInput): Complex { const zc = Complex.from(z); return new Complex(this.re + zc.re, this.im + zc.im); }
  sub(z: ComplexInput): Complex { const zc = Complex.from(z); return new Complex(this.re - zc.re, this.im - zc.im); }
  neg(): Complex { return new Complex(-this.re, -this.im); }

  mul(z: ComplexInput): Complex {
    const zc = Complex.from(z);
    return new Complex(this.re * zc.re - this.im * zc.im, this.re * zc.im + this.im * zc.re);
  }

  div(z: ComplexInput): Complex {
    const zc = Complex.from(z);
    const d = zc.re * zc.re + zc.im * zc.im;
    if (d === 0) throw new RangeError('Complex division by zero');
    return new Complex((this.re * zc.re + this.im * zc.im) / d, (this.im * zc.re - this.re * zc.im) / d);
  }

  inv(): Complex {
    const d = this.re * this.re + this.im * this.im;
    if (d === 0) throw new RangeError('Complex inversion of zero');
    return new Complex(this.re / d, -this.im / d);
  }

  conj(): Complex { return new Complex(this.re, -this.im); }
  abs2(): number { return this.re * this.re + this.im * this.im; }
  abs(): number { return Math.hypot(this.re, this.im); }
  arg(): number { return Math.atan2(this.im, this.re); }

  powReal(t: number): Complex {
    const r = this.abs();
    const theta = this.arg();
    const rt = Math.pow(r, t);
    return new Complex(rt * Math.cos(theta * t), rt * Math.sin(theta * t));
  }

  equals(z: ComplexInput, eps = EPS): boolean {
    const zc = Complex.from(z);
    return this.sub(zc).abs() <= eps;
  }

  clean(eps = EPS): Complex {
    const re = Math.abs(this.re) < eps ? 0 : this.re;
    const im = Math.abs(this.im) < eps ? 0 : this.im;
    return new Complex(re, im);
  }

  toString(precision = 15): string {
    const z = this.clean();
    const re = Number(z.re.toPrecision(precision));
    const im = Number(z.im.toPrecision(precision));
    if (im === 0) return String(re);
    if (re === 0) return `${im}im`;
    return `${re}${im < 0 ? '' : '+'}${im}im`;
  }
}

export const C = (re = 0, im = 0): Complex => new Complex(re, im);
export const I = new Complex(0, 1);

export function cexpi(theta: number): Complex {
  return new Complex(Math.cos(theta), Math.sin(theta));
}

export function cdist(a: ComplexInput, b: ComplexInput): number {
  return Complex.from(a).sub(b).abs();
}
