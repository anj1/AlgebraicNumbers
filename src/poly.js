import { allRoots, allRootsCertifiedSimplified } from 'flo-poly';
import { Complex, C, EPS } from './complex.js';

export function trimAscending(p, eps = EPS) {
  const q = p.map(Number);
  while (q.length > 1 && Math.abs(q[q.length - 1]) <= eps) q.pop();
  return q.length ? q : [0];
}

export function trimDescending(p, eps = EPS) {
  const q = p.map(Number);
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
  const lc = p[p.length - 1];
  if (Math.abs(lc) <= eps) return [0];
  return cleanRealCoeffs(p.map(c => c / lc), eps);
}

export function cleanRealCoeffs(p, eps = 1e-9) {
  return trimAscending(p.map(c => {
    if (Math.abs(c) < eps) return 0;
    const r = Math.round(c);
    return Math.abs(c - r) < eps ? r : c;
  }), eps);
}

export function addPoly(p, q) {
  const n = Math.max(p.length, q.length);
  const r = Array(n).fill(0);
  for (let i = 0; i < n; i++) r[i] = (p[i] ?? 0) + (q[i] ?? 0);
  return cleanRealCoeffs(r);
}

export function subPoly(p, q) {
  const n = Math.max(p.length, q.length);
  const r = Array(n).fill(0);
  for (let i = 0; i < n; i++) r[i] = (p[i] ?? 0) - (q[i] ?? 0);
  return cleanRealCoeffs(r);
}

export function scalePoly(p, a) {
  return cleanRealCoeffs(p.map(c => a * c));
}

export function mulPoly(p, q) {
  const r = Array(p.length + q.length - 1).fill(0);
  for (let i = 0; i < p.length; i++) {
    for (let j = 0; j < q.length; j++) r[i + j] += p[i] * q[j];
  }
  return cleanRealCoeffs(r);
}

export function divPolyExact(numer, denom, eps = 1e-8) {
  numer = trimAscending(numer, eps).slice();
  denom = trimAscending(denom, eps);
  if (degree(denom) < 0 || Math.abs(denom[denom.length - 1]) <= eps) throw new RangeError('division by zero polynomial');
  if (degree(numer) < degree(denom)) return { quotient: [0], remainder: numer };
  const q = Array(degree(numer) - degree(denom) + 1).fill(0);
  const ddeg = degree(denom);
  const dlc = denom[denom.length - 1];
  for (let k = degree(numer) - ddeg; k >= 0; k--) {
    const coeff = numer[ddeg + k] / dlc;
    q[k] = coeff;
    for (let j = 0; j <= ddeg; j++) numer[j + k] -= coeff * denom[j];
  }
  const rem = cleanRealCoeffs(numer.slice(0, ddeg), eps);
  return { quotient: cleanRealCoeffs(q, eps), remainder: rem };
}

export function derivativeAscending(p) {
  if (p.length <= 1) return [0];
  const r = [];
  for (let i = 1; i < p.length; i++) r.push(i * p[i]);
  return cleanRealCoeffs(r);
}

export function evalPolyAscending(p, z) {
  z = Complex.from(z);
  let acc = C(0, 0);
  for (let i = p.length - 1; i >= 0; i--) acc = acc.mul(z).add(p[i]);
  return acc;
}

export function evalPolyRealAscending(p, x) {
  let acc = 0;
  for (let i = p.length - 1; i >= 0; i--) acc = acc * x + p[i];
  return acc;
}

export function allRealRootsAscending(p) {
  p = trimAscending(p);
  if (degree(p) <= 0) return [];
  const roots = allRoots(ascendingToDescending(p));
  return roots.map(rootToNumber).filter(Number.isFinite).sort((a, b) => a - b);
}

export function allRealRootsCertifiedAscending(p) {
  p = trimAscending(p);
  if (degree(p) <= 0) return [];
  return allRootsCertifiedSimplified(ascendingToDescending(p));
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
  if (n === 1) return [C(-p[0] / p[1], 0)];

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
  if (n === 1) return [C(-p[0] / p[1], 0)];

  const maxIter = opts.maxIter ?? 2000;
  const tol = opts.tol ?? 1e-12;
  const radius = 1 + Math.max(...p.slice(0, -1).map(Math.abs));
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
  return polyFromRoots(roots);
}

export function composedProduct(p, q) {
  const rp = rootsAscending(p);
  const rq = rootsAscending(q);
  const roots = [];
  for (const a of rp) for (const b of rq) roots.push(a.mul(b));
  return polyFromRoots(roots);
}

export function interleaveZerosAscending(p, nZeros) {
  if (nZeros < 0) throw new RangeError('nZeros must be non-negative');
  const out = [];
  for (let i = 0; i < p.length; i++) {
    out.push(p[i]);
    if (i !== p.length - 1) for (let j = 0; j < nZeros; j++) out.push(0);
  }
  return trimAscending(out);
}

export function syntheticDivideByLinear(p, root, eps = 1e-8) {
  p = trimAscending(p, eps);
  const n = degree(p);
  if (n <= 0) return { quotient: [0], remainder: p[0] ?? 0 };
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
  return p.every((c, i) => Math.abs(c - q[i]) <= eps * Math.max(1, Math.abs(c), Math.abs(q[i])));
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
