import test from 'node:test';
import assert from 'node:assert/strict';
import {
  AlgebraicNumber,
  allRealRootsAscending,
  algRoots,
  confirmAlgebraicNumber,
  cosAlg,
  factorPolynomial,
  invTotient,
  root,
  samePolynomialUpToScale,
  sqrt
} from '../src/index.js';

const near = (a, b, eps = 1e-7) => assert.ok(Math.abs(a - b) < eps, `${a} not near ${b}`);
const nearComplex = (z, re, im = 0, eps = 1e-7) => {
  near(z.re, re, eps);
  near(z.im, im, eps);
};

test('flo-poly real roots are wired with ascending coefficients', () => {
  const p = [720, -1764, 1624, -735, 175, -21, 1];
  const roots = allRealRootsAscending(p);
  assert.equal(roots.length, 6);
  roots.forEach((r, i) => near(r, i + 1, 1e-6));
});

test('constructs rational/integer algebraic numbers', () => {
  const a = AlgebraicNumber.from(7);
  assert.deepEqual(a.coeff, [-7, 1]);
  nearComplex(a.apprx, 7, 0);

  const b = AlgebraicNumber.rational(3, 4);
  nearComplex(b.apprx, 0.75, 0);
  assert.ok(samePolynomialUpToScale(b.coeff, [-0.75, 1]));
});

test('adds, multiplies, divides', () => {
  const two = AlgebraicNumber.from(2);
  const three = AlgebraicNumber.from(3);

  nearComplex(two.add(three).apprx, 5);
  nearComplex(three.mul(two).apprx, 6);
  nearComplex(three.div(two).apprx, 1.5);

  assert.ok(two.add(three).equals(5));
  assert.ok(three.mul(two).equals(6));
});

test('radicals preserve selected root and simplify reducible cases', () => {
  const s2 = sqrt(AlgebraicNumber.from(2));
  nearComplex(s2.apprx, Math.SQRT2);
  assert.ok(samePolynomialUpToScale(s2.coeff, [-2, 0, 1]));

  const s4 = sqrt(AlgebraicNumber.from(4));
  nearComplex(s4.apprx, 2);
  assert.ok(samePolynomialUpToScale(s4.coeff, [-2, 1]));
});

test('complex algebraic number i squares to -1', () => {
  const i = root(AlgebraicNumber.from(-1), 2);
  nearComplex(i.apprx, 0, 1);
  assert.ok(samePolynomialUpToScale(i.coeff, [1, 0, 1]));

  const minusOne = i.mul(i);
  nearComplex(minusOne.apprx, -1, 0);
  assert.ok(minusOne.equals(-1));
});

test('algRoots returns algebraic roots of a polynomial', () => {
  const roots = algRoots([-6, 11, -6, 1]).map(x => x.apprx.re).sort((a, b) => a - b);
  roots.forEach((r, i) => near(r, i + 1));
});

test('algebraic-numbers-backed factorization falls back to rational-root factors', () => {
  const factors = factorPolynomial([-4, 0, 1]);
  assert.equal(factors.length, 2);
  assert.ok(factors.some(f => samePolynomialUpToScale(f, [-2, 1])));
  assert.ok(factors.some(f => samePolynomialUpToScale(f, [2, 1])));
});

test('trigonometric algebraic helpers', () => {
  const c = cosAlg(1, 3);
  nearComplex(c.apprx, 0.5, 0, 1e-7);
  assert.ok(c.equals(AlgebraicNumber.rational(1, 2)));
});

test('inverse totient helper', () => {
  assert.deepEqual(invTotient(2), [3, 4, 6]);
});

test('residual is near zero', () => {
  const s3 = sqrt(AlgebraicNumber.from(3));
  const r = confirmAlgebraicNumber(s3);
  nearComplex(r, 0, 0, 1e-7);
});
