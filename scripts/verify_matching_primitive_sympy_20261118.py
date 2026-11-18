"""Independent symbolic identities for two alternative cold-EOS interpolants."""
import json
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
out=ROOT/"data/matching_interpolation_20261118/SYMBOLIC_VERIFICATION.json"
t,q,a,d,alpha=sp.symbols("t q a d alpha",positive=True)
m=1-q
phi=t/(q+m*t)
F0=t/m-q*sp.log(1+m*t/q)/m**2
F1=t**2/(2*m)-q*t/m**2+q*q*sp.log(1+m*t/q)/m**3
checks={
 "rational_F0_derivative":sp.simplify(sp.diff(F0,t)-phi),
 "rational_F1_derivative":sp.simplify(sp.diff(F1,t)-t*phi),
 "power_primitive":sp.simplify(sp.diff(a*t**(alpha+1)/(alpha+1)+d*t**(alpha+2)/(alpha+2),t)-(a+d*t)*t**alpha),
}
receipt=dict(tool="system Python SymPy",version=sp.__version__,identities={k:str(v) for k,v in checks.items()},all_zero=all(v==0 for v in checks.values()),
  caveat="Rational expression requires q>0,q!=1; code uses continuous q=1 treatment. Does not establish physical superiority or whole-EOS convergence.")
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
raise SystemExit(0 if receipt["all_zero"] else 1)
