export function primeFactorization(n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError('n must be a positive integer');
  const out = new Map();
  let m = n;
  for (let p = 2; p * p <= m; p += p === 2 ? 1 : 2) {
    while (m % p === 0) {
      out.set(p, (out.get(p) ?? 0) + 1);
      m = Math.trunc(m / p);
    }
  }
  if (m > 1) out.set(m, (out.get(m) ?? 0) + 1);
  return out;
}

export function allDivisorsFromFactors(factors) {
  const entries = [...factors.entries()];
  const out = [];

  function rec(i, cur) {
    if (i === entries.length) { out.push(cur); return; }
    const [p, e] = entries[i];
    let m = 1;
    for (let k = 0; k <= e; k++) {
      rec(i + 1, cur * m);
      m *= p;
    }
  }

  rec(0, 1);
  return out.sort((a, b) => a - b);
}

export function divisors(n) {
  return allDivisorsFromFactors(primeFactorization(n));
}

export function isPrime(n) {
  if (!Number.isInteger(n) || n < 2) return false;
  if (n === 2) return true;
  if (n % 2 === 0) return false;
  for (let p = 3; p * p <= n; p += 2) if (n % p === 0) return false;
  return true;
}

export function eulerPhi(n) {
  let out = n;
  for (const p of primeFactorization(n).keys()) out = Math.trunc((out / p) * (p - 1));
  return out;
}

export function invTotient(x) {
  if (!Number.isInteger(x) || x < 1) throw new RangeError('x must be a positive integer');
  const invs = new Set();

  function totientReps(rem, factorList) {
    if (rem < 1) return;
    if (rem === 1) {
      let m = 1;
      for (let i = 1; i < factorList.length; i++) {
        const [p, gamma] = factorList[i];
        m *= p ** (gamma + 1);
      }
      invs.add(m);
    }

    const factors = primeFactorization(rem);

    for (const divisor of allDivisorsFromFactors(factors)) {
      const d = divisor + 1;
      if (isPrime(d) && d > factorList[factorList.length - 1][0]) {
        totientReps(Math.trunc(rem / divisor), [...factorList, [d, 0]]);
      }
    }

    for (const [p, e] of factors.entries()) {
      for (let gamma = 1; gamma <= e; gamma++) {
        const d = (p ** gamma) * (p - 1);
        if (rem % d === 0 && p > factorList[factorList.length - 1][0]) {
          totientReps(Math.trunc(rem / d), [...factorList, [p, gamma]]);
        }
      }
    }
  }

  totientReps(x, [[0, 0]]);
  return [...invs].sort((a, b) => a - b);
}
