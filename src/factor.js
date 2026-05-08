import algebrite from 'algebrite';
import {
  addPoly,
  cleanRealCoeffs,
  degree,
  divPolyExact,
  mulPoly,
  normalizeAscending,
  samePolynomialUpToScale,
  syntheticDivideByLinear,
  trimAscending,
  allRealRootsAscending
} from './poly.js';

const EPS = 1e-8;

export function polynomialToalgebriteExpr(p, variable = 'x') {
  p = trimAscending(p);
  if (p.length === 1) return String(p[0]);
  const terms = [];
  for (let i = p.length - 1; i >= 0; i--) {
    const c = p[i];
    if (Math.abs(c) < EPS) continue;
    const sign = c < 0 ? '-' : '+';
    const a = Math.abs(c);
    let body;
    if (i === 0) body = formatNumber(a);
    else if (i === 1) body = a === 1 ? variable : `${formatNumber(a)}*${variable}`;
    else body = a === 1 ? `${variable}^${i}` : `${formatNumber(a)}*${variable}^${i}`;
    terms.push({ sign, body });
  }
  if (terms.length === 0) return '0';
  return terms.map((t, i) => i === 0 && t.sign === '+' ? t.body : `${t.sign}${t.body}`).join('');
}

function formatNumber(x) {
  const r = Math.round(x);
  if (Math.abs(x - r) < EPS) return String(r);
  return Number(x.toPrecision(15)).toString();
}

export function factorWithalgebrite(p, variable = 'x') {
  const expr = polynomialToalgebriteExpr(p, variable);
  const factored = String(algebrite.run(`factor(${expr})`));
  const parsed = parseFactoredPolynomial(factored, variable);
  if (parsed.length <= 1) return [];

  const product = parsed.reduce((acc, f) => mulPoly(acc, f), [1]);
  return samePolynomialUpToScale(product, p, 1e-6) ? parsed.map(f => normalizeAscending(f)) : [];
}

export function factorPolynomial(p) {
  p = normalizeAscending(p);
  if (degree(p) <= 1) return [p];

  try {
    const factors = factorWithalgebrite(p);
    if (factors.length > 1) return factors;
  } catch {
    // Fall through to a numeric rational-root factorization. algebrite output formatting
    // can vary, but it is still the primary symbolic factorization path.
  }

  return factorByRationalRoots(p);
}

export function parseFactoredPolynomial(src, variable = 'x') {
  const s = insertImplicitMultiplication(String(src).replace(/\s+/g, ''));
  const rawFactors = splitTopLevelProduct(s);
  if (rawFactors.length <= 1) return [];
  const out = [];
  for (const raw of rawFactors) {
    const { base, exponent } = peelPower(raw);
    const poly = parsePolynomialExpression(stripOuterParens(base), variable);
    if (!poly) return [];
    for (let i = 0; i < exponent; i++) out.push(cleanRealCoeffs(poly));
  }
  return out;
}

function insertImplicitMultiplication(s) {
  return s
    .replace(/\)\(/g, ')*(')
    .replace(/(\d|x)\(/g, '$1*(')
    .replace(/\)(x|\d)/g, ')*$1')
    .replace(/(\d)x/g, '$1*x');
}

function splitTopLevelProduct(s) {
  const out = [];
  let depth = 0;
  let start = 0;
  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    if (ch === '(') depth++;
    else if (ch === ')') depth--;
    else if (ch === '*' && depth === 0) {
      out.push(s.slice(start, i));
      start = i + 1;
    }
  }
  out.push(s.slice(start));
  return out.filter(Boolean);
}

function peelPower(raw) {
  let s = raw;
  let exponent = 1;
  const m = s.match(/^(.*)\^(\d+)$/);
  if (m) {
    s = m[1];
    exponent = Number(m[2]);
  }
  return { base: s, exponent };
}

function stripOuterParens(s) {
  while (s.startsWith('(') && s.endsWith(')')) {
    let depth = 0;
    let wraps = true;
    for (let i = 0; i < s.length; i++) {
      if (s[i] === '(') depth++;
      if (s[i] === ')') depth--;
      if (depth === 0 && i < s.length - 1) { wraps = false; break; }
    }
    if (!wraps) break;
    s = s.slice(1, -1);
  }
  return s;
}

function parsePolynomialExpression(input, variable = 'x') {
  const tokens = tokenize(input, variable);
  let pos = 0;

  function peek() { return tokens[pos]; }
  function take(t) { if (peek() === t) { pos++; return true; } return false; }

  function parseExpr() {
    let acc = parseTerm();
    while (peek() === '+' || peek() === '-') {
      const op = tokens[pos++];
      const rhs = parseTerm();
      acc = op === '+' ? addPoly(acc, rhs) : addPoly(acc, rhs.map(c => -c));
    }
    return acc;
  }

  function parseTerm() {
    let acc = parsePower();
    while (take('*')) acc = mulPoly(acc, parsePower());
    return acc;
  }

  function parsePower() {
    let base = parseUnary();
    if (take('^')) {
      const e = tokens[pos++];
      if (typeof e !== 'number' || !Number.isInteger(e) || e < 0) throw new Error('polynomial exponents must be non-negative integers');
      let acc = [1];
      for (let i = 0; i < e; i++) acc = mulPoly(acc, base);
      base = acc;
    }
    return base;
  }

  function parseUnary() {
    if (take('+')) return parseUnary();
    if (take('-')) return parseUnary().map(c => -c);
    return parsePrimary();
  }

  function parsePrimary() {
    const t = tokens[pos++];
    if (typeof t === 'number') return [t];
    if (t === variable) return [0, 1];
    if (t === '(') {
      const e = parseExpr();
      if (!take(')')) throw new Error('missing )');
      return e;
    }
    throw new Error(`unexpected token ${t}`);
  }

  const poly = cleanRealCoeffs(parseExpr());
  if (pos !== tokens.length) throw new Error(`unparsed tokens: ${tokens.slice(pos).join(' ')}`);
  return poly;
}

function tokenize(s, variable) {
  const tokens = [];
  let i = 0;
  while (i < s.length) {
    const ch = s[i];
    if ('+-*^()'.includes(ch)) { tokens.push(ch); i++; continue; }
    if (ch === variable) { tokens.push(variable); i++; continue; }
    if (/\d|\./.test(ch)) {
      let j = i + 1;
      while (j < s.length && /[\d.eE+-]/.test(s[j])) {
        if ((s[j] === '+' || s[j] === '-') && !/[eE]/.test(s[j - 1])) break;
        j++;
      }
      const num = Number(s.slice(i, j));
      if (!Number.isFinite(num)) throw new Error(`bad number ${s.slice(i, j)}`);
      tokens.push(num);
      i = j;
      continue;
    }
    throw new Error(`bad char ${ch}`);
  }
  return tokens;
}

function rationalApprox(x, maxDen = 256, eps = 1e-8) {
  let best = null;
  for (let q = 1; q <= maxDen; q++) {
    const p = Math.round(x * q);
    const err = Math.abs(x - p / q);
    if (err < eps && (!best || q < best.den)) best = { num: p, den: q, value: p / q };
  }
  return best;
}

function factorByRationalRoots(p) {
  let residual = normalizeAscending(p);
  const factors = [];
  let changed = true;

  while (changed && degree(residual) > 1) {
    changed = false;
    let candidates = [];
    try { candidates = [...new Set(allRealRootsAscending(residual).map(r => rationalApprox(r)).filter(Boolean).map(r => r.value))]; }
    catch { candidates = []; }

    // Deterministic small-rational backup for cases where real-root isolation misses a
    // repeated root after numerical deflation.
    for (let q = 1; q <= 32; q++) {
      for (let n = -128; n <= 128; n++) candidates.push(n / q);
    }

    for (const r of candidates) {
      const { quotient, remainder } = syntheticDivideByLinear(residual, r);
      if (Math.abs(remainder) < 1e-6) {
        factors.push(cleanRealCoeffs([-r, 1]));
        residual = normalizeAscending(quotient);
        changed = true;
        break;
      }
    }
  }

  if (degree(residual) >= 1) factors.push(residual);
  return factors.length ? factors : [normalizeAscending(p)];
}
