/*eslint no-constant-condition: ["error", { "checkLoops": false }]*/
"use strict";
function normcdf(z, mean = 0.0, sd = 1.0) {
  // implementation of ZKERCOOK algorithm from Brophy and Wood, 1989.
  const x = (z - mean) / sd;
  let i = 2, p = 1.0, p1 = 0.0, t = 1.0, t1 = x*x/4.0, z1 = t1;
  if(Math.abs(x) > 6.8) {
    const zz = x * x;
    const w  = 1.0 - 1.0 / (zz + 3.0 - 1. / (0.22 * zz + 0.704));
    p1 = 1 / Math.sqrt(2.0 * Math.PI) * Math.exp(-zz / 2.0) * w / x;
    if(Math.sign(x) === -1) { p1 = 1 - p1; }
  } else {
    while(true)
    {
      t = z1 * (t1 - t) / i;
      i++;
      p = p + t / i;

      if(p == p1) { break; }
      
      t1 = z1 * (t - t1) / i;
      i++;
      p1 = p;
    }
    p = 1 / Math.sqrt(2.0 * Math.PI) * x * Math.exp(-z1 / 2.0) * p;
    p1 = 0.5 - p;
  }
  return 1.0 - p1
}

function erf(z) {
  // https://stats.stackexchange.com/a/187909/311838
  return 2.0 * normcdf(z * Math.sqrt(2)) - 1;
}

function normpdf(x, mean = 0.0, sd = 1.0) {
  return Math.exp(-1.0/2.0 * ((x - mean)/sd)**2) / (sd * Math.sqrt(2 * Math.PI))
}

function tcdf(x, v) {
  // Li and De Moor (1999), by way of Yerukala, Boiroju, and Reddy (2013)
  if(!Number.isInteger(v) || v < 1) { return null; }
  if(v >= 30) { return normcdf(x); }
  if(v == 1) { return 0.5 + 1.0 / Math.PI * Math.atan(x); }
  if(v == 2) { return 0.5 + x / (2.0 * Math.sqrt(2+x*x)); }
  const l = (4*v + x*x - 1.0)/(4.0*v + 2.0*x*x);
  return normcdf(x * l);
}
