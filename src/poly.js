import { allRoots, allRootsCertifiedSimplified } from 'flo-poly';
import { Complex, C, EPS } from './complex.js';

function gcd(a, b) {
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

export function trimAscending(p, eps = EPS) {
  const q = p.slice();
  const isBig = typeof q[0] === 'bigint';
  if (isBig) {
    while (q.length > 1 && q[q.length - 1] === 0n) q.pop();
    return q.length ? q : [0n];
  }
  while (q.length > 1 && Math.abs(q[q.length - 1]) <= eps) q.pop();
  return q.length ? q : [0];
}

export function trimDescending(p, eps = EPS) {
  const q = p.slice();
  const isBig = typeof q[0] === 'bigint';
  if (isBig) {
    while (q.length > 1 && q[0] === 0n) q.shift();
    return q.length ? q : [0n];
  }
  while (q.length > 1 && Math.abs(q[0]) <= eps) q.shift();
  return q.length ? q : [0];
}

export function ascendingToDescending(p) {
  return trimAscending(p).slice().reverse();
}

export function descendingToAscending(p) {
  return trimDescending(p).slice().reverse();
}

export function degree(p) {
  return trimAscending(p).length - 1;
}

export function normalizeAscending(p, eps = EPS) {
  p = trimAscending(p, eps);
  if (p.length === 0) return [typeof p[0] === 'bigint' ? 0n : 0];
  const lc = p[p.length - 1];
  if (typeof lc === 'bigint') {
    if (lc === 0n) return [0n];
    let g = p[0] < 0n ? -p[0] : p[0];
    for (let i = 1; i < p.length; i++) {
       const c = p[i] < 0n ? -p[i] : p[i];
       g = g === 0n ? c : gcd(g, c);
    }
    if (g === 0n) return [0n];
    const sign = lc < 0n ? -1n : 1n;
    return p.map(c => (c * sign) / g);
  }
  if (Math.abs(lc) <= eps) return [0];
  return cleanRealCoeffs(p.map(c => c / lc), eps);
}

export function cleanRealCoeffs(p, eps = 1e-9) {
  if (p.length > 0 && typeof p[0] === 'bigint') return trimAscending(p);
  return trimAscending(p.map(c => {
    if (Math.abs(c) < eps) return 0;
    const r = Math.round(c);
    return Math.abs(c - r) < eps ? r : c;
  }), eps);
}

export function addPoly(p, q) {
  const n = Math.max(p.length, q.length);
  const isBig = typeof (p[0] ?? q[0]) === 'bigint';
  const r = Array(n).fill(isBig ? 0n : 0);
  for (let i = 0; i < n; i++) {
     let a = p[i] ?? (isBig ? 0n : 0);
     let b = q[i] ?? (isBig ? 0n : 0);
     if (isBig) { a = BigInt(a); b = BigInt(b); }
     r[i] = a + b;
  }
  return isBig ? trimAscending(r) : cleanRealCoeffs(r);
}

export function subPoly(p, q) {
  const n = Math.max(p.length, q.length);
  const isBig = typeof (p[0] ?? q[0]) === 'bigint';
  const r = Array(n).fill(isBig ? 0n : 0);
  for (let i = 0; i < n; i++) {
     let a = p[i] ?? (isBig ? 0n : 0);
     let b = q[i] ?? (isBig ? 0n : 0);
     if (isBig) { a = BigInt(a); b = BigInt(b); }
     r[i] = a - b;
  }
  return isBig ? trimAscending(r) : cleanRealCoeffs(r);
}

export function scalePoly(p, a) {
  const isBig = typeof p[0] === 'bigint';
  if (isBig) return trimAscending(p.map(c => BigInt(Math.round(Number(a))) * c));
  return cleanRealCoeffs(p.map(c => a * c));
}

export function mulPoly(p, q) {
  if (p.length === 0 || q.length === 0) return [];
  const isBig = typeof (p[0] ?? q[0]) === 'bigint';
  const r = Array(p.length + q.length - 1).fill(isBig ? 0n : 0);
  for (let i = 0; i < p.length; i++) {
    for (let j = 0; j < q.length; j++) {
       if (isBig) r[i + j] += BigInt(p[i]) * BigInt(q[j]);
       else r[i + j] += p[i] * q[j];
    }
  }
  return isBig ? trimAscending(r) : cleanRealCoeffs(r);
}

export function divPolyExact(numer, denom, eps = 1e-8) {
  numer = trimAscending(numer, eps).slice();
  denom = trimAscending(denom, eps);
  const isBig = typeof numer[0] === 'bigint' || typeof denom[0] === 'bigint';
  if (degree(denom) < 0 || (isBig ? denom[denom.length - 1] === 0n : Math.abs(denom[denom.length - 1]) <= eps)) throw new RangeError('division by zero polynomial');
  if (degree(numer) < degree(denom)) return { quotient: [isBig ? 0n : 0], remainder: numer };
  const q = Array(degree(numer) - degree(denom) + 1).fill(isBig ? 0n : 0);
  const ddeg = degree(denom);
  const dlc = denom[denom.length - 1];
  for (let k = degree(numer) - ddeg; k >= 0; k--) {
    if (isBig && numer[ddeg + k] % dlc !== 0n) throw new Error('not exactly divisible');
    const coeff = isBig ? numer[ddeg + k] / dlc : numer[ddeg + k] / dlc;
    q[k] = coeff;
    for (let j = 0; j <= ddeg; j++) numer[j + k] -= coeff * denom[j];
  }
  const rem = isBig ? trimAscending(numer.slice(0, ddeg)) : cleanRealCoeffs(numer.slice(0, ddeg), eps);
  return { quotient: isBig ? trimAscending(q) : cleanRealCoeffs(q, eps), remainder: rem };
}

export function derivativeAscending(p) {
  if (p.length <= 1) return typeof p[0] === 'bigint' ? [0n] : [0];
  const isBig = typeof p[0] === 'bigint';
  const r = [];
  for (let i = 1; i < p.length; i++) r.push(isBig ? BigInt(i) * p[i] : i * p[i]);
  return isBig ? trimAscending(r) : cleanRealCoeffs(r);
}

export function evalPolyAscending(p, z) {
  z = Complex.from(z);
  let acc = C(0, 0);
  for (let i = p.length - 1; i >= 0; i--) acc = acc.mul(z).add(Number(p[i]));
  return acc;
}

export function evalPolyRealAscending(p, x) {
  let acc = 0;
  for (let i = p.length - 1; i >= 0; i--) acc = acc * x + Number(p[i]);
  return acc;
}

export function allRealRootsAscending(p) {
  p = trimAscending(p);
  if (degree(p) <= 0) return [];
  const roots = allRoots(ascendingToDescending(p).map(Number));
  return roots.map(rootToNumber).filter(Number.isFinite).sort((a, b) => a - b);
}

export function allRealRootsCertifiedAscending(p) {
  p = trimAscending(p);
  if (degree(p) <= 0) return [];
  return allRootsCertifiedSimplified(ascendingToDescending(p).map(Number));
}

function rootToNumber(r) {
  if (typeof r === 'number') return r;
  if (Array.isArray(r) && r.length >= 2) return (Number(r[0]) + Number(r[1])) / 2;
  if (r && typeof r.tS === 'number' && typeof r.tE === 'number') return (r.tS + r.tE) / 2;
  if (r && typeof r.mid === 'number') return r.mid;
  return Number(r);
}

export function rootsAscending(p, opts = {}) {
  p = normalizeAscending(p);
  const n = degree(p);
  if (n <= 0) return [];
  if (n === 1) return [C(-Number(p[0]) / Number(p[1]), 0)];

  try {
    const real = allRealRootsAscending(p);
    if (real.length === n) return real.map(x => C(x, 0));
  } catch {
    // flo-poly intentionally only targets real roots; use the complex fallback below.
  }

  return durandKernerRoots(p, opts);
}

export function durandKernerRoots(p, opts = {}) {
  p = normalizeAscending(p);
  const n = degree(p);
  if (n <= 0) return [];
  if (n === 1) return [C(-Number(p[0]) / Number(p[1]), 0)];

  const maxIter = opts.maxIter ?? 2000;
  const tol = opts.tol ?? 1e-12;
  const radius = 1 + Math.max(...p.slice(0, -1).map(c => Math.abs(Number(c))));
  let roots = Array.from({ length: n }, (_, k) => {
    const theta = (2 * Math.PI * k) / n + 0.2718281828459045;
    return C(radius * Math.cos(theta), radius * Math.sin(theta));
  });

  for (let iter = 0; iter < maxIter; iter++) {
    let maxDelta = 0;
    const next = roots.map((zk, k) => {
      let denom = C(1, 0);
      for (let j = 0; j < n; j++) if (j !== k) denom = denom.mul(zk.sub(roots[j]));
      if (denom.abs() === 0) denom = C(tol, tol);
      const delta = evalPolyAscending(p, zk).div(denom);
      maxDelta = Math.max(maxDelta, delta.abs());
      return zk.sub(delta);
    });
    roots = next;
    if (maxDelta < tol) break;
  }

  return roots.map(z => z.clean(1e-9)).sort((a, b) => a.re === b.re ? a.im - b.im : a.re - b.re);
}

export function polyFromRoots(roots, eps = 1e-8) {
  let coeffs = [C(1, 0)];
  for (const root of roots.map(Complex.from)) {
    const next = Array(coeffs.length + 1).fill(null).map(() => C(0, 0));
    for (let i = 0; i < coeffs.length; i++) {
      next[i] = next[i].sub(coeffs[i].mul(root));
      next[i + 1] = next[i + 1].add(coeffs[i]);
    }
    coeffs = next;
  }
  const real = coeffs.map(c => {
    if (Math.abs(c.im) > eps) throw new Error(`non-real coefficient ${c.toString()}`);
    return c.re;
  });
  return normalizeAscending(cleanRealCoeffs(real, eps));
}

export function composedSum(p, q) {
  const rp = rootsAscending(p);
  const rq = rootsAscending(q);
  const roots = [];
  for (const a of rp) for (const b of rq) roots.push(a.add(b));
  const f = polyFromRoots(roots);
  const isBig = typeof p[0] === 'bigint' || typeof q[0] === 'bigint';
  if (isBig) {
     let lc1 = p[p.length - 1]; if (typeof lc1 !== 'bigint') lc1 = BigInt(lc1);
     let lc2 = q[q.length - 1]; if (typeof lc2 !== 'bigint') lc2 = BigInt(lc2);
     const lc = (lc1 ** BigInt(degree(q))) * (lc2 ** BigInt(degree(p)));
     const numLc = Number(lc);
     return normalizeAscending(f.map(c => BigInt(Math.round(c * numLc))));
  }
  return f;
}

export function composedProduct(p, q) {
  const rp = rootsAscending(p);
  const rq = rootsAscending(q);
  const roots = [];
  for (const a of rp) for (const b of rq) roots.push(a.mul(b));
  const f = polyFromRoots(roots);
  const isBig = typeof p[0] === 'bigint' || typeof q[0] === 'bigint';
  if (isBig) {
     let lc1 = p[p.length - 1]; if (typeof lc1 !== 'bigint') lc1 = BigInt(lc1);
     let lc2 = q[q.length - 1]; if (typeof lc2 !== 'bigint') lc2 = BigInt(lc2);
     const lc = (lc1 ** BigInt(degree(q))) * (lc2 ** BigInt(degree(p)));
     const numLc = Number(lc);
     return normalizeAscending(f.map(c => BigInt(Math.round(c * numLc))));
  }
  return f;
}

export function interleaveZerosAscending(p, nZeros) {
  if (nZeros < 0) throw new RangeError('nZeros must be non-negative');
  const out = [];
  const isBig = typeof p[0] === 'bigint';
  const zero = isBig ? 0n : 0;
  for (let i = 0; i < p.length; i++) {
    out.push(p[i]);
    if (i !== p.length - 1) for (let j = 0; j < nZeros; j++) out.push(zero);
  }
  return trimAscending(out);
}

export function syntheticDivideByLinear(p, root, eps = 1e-8) {
  p = trimAscending(p, eps);
  const n = degree(p);
  const isBig = typeof p[0] === 'bigint';
  if (n <= 0) return { quotient: [isBig ? 0n : 0], remainder: p[0] ?? (isBig ? 0n : 0) };
  
  if (isBig) p = p.map(Number);
  const q = Array(n).fill(0);
  q[n - 1] = p[n];
  for (let k = n - 1; k >= 1; k--) q[k - 1] = p[k] + root * q[k];
  const rem = p[0] + root * q[0];
  return { quotient: cleanRealCoeffs(q, eps), remainder: rem };
}

export function samePolynomialUpToScale(p, q, eps = 1e-8) {
  p = normalizeAscending(p, eps);
  q = normalizeAscending(q, eps);
  if (p.length !== q.length) return false;
  const isBig = typeof p[0] === 'bigint' && typeof q[0] === 'bigint';
  if (isBig) return p.every((c, i) => c === q[i]);
  return p.every((c, i) => Math.abs(Number(c) - Number(q[i])) <= eps * Math.max(1, Math.abs(Number(c)), Math.abs(Number(q[i]))));
}

export function minRootSeparation(p) {
  const roots = rootsAscending(p);
  if (roots.length <= 1) return Infinity;
  let min = Infinity;
  for (let i = 0; i < roots.length; i++) {
    for (let j = i + 1; j < roots.length; j++) min = Math.min(min, roots[i].sub(roots[j]).abs());
  }
  return min;
}
