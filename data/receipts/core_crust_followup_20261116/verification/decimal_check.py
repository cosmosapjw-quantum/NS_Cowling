"""Read-only independent Decimal oracle of spherical interface helpers.
Run: PYTHONPATH=python_frontend python3 build/review_scurto_peer/check.py
Tests exact represented inputs, not parameter-fit or CLD admissibility.
"""
import hashlib, json, math, pathlib, platform, subprocess
from decimal import Decimal as D, localcontext
from compose_eos.scurto_interface import spherical_phi, surface_sigma, surface_sigma_prime

ys = sorted(set([1e-8, 1e-7, 1e-5, .01, .1, .3, .49, .49999999,
     math.nextafter(.5, 0), .5, math.nextafter(.5, 1), .50000001,
     .51, .7, .9, .99, 1-1e-5, 1-1e-7, 1-1e-8]))
bs = [math.nextafter(-16, math.inf), -15.999999, -15., -1., 0., 7., 1e10, 1e50, 1e100]
ss = [0., math.nextafter(0, 1), 1e-300, 1.3, 1e6]
fs = sorted(set([0., math.nextafter(0, 1), 1e-300, 1e-30, 1e-8, .001, .125, .5, .9, 1-1e-8, math.nextafter(1, 0), 1.]))
rows=[]
with localcontext() as ctx:
    ctx.prec=400
    for y in ys:
        for b in bs:
            for s in ss:
                x, bb, ss0 = map(D.from_float,(y,b,s))
                dd=bb+x**-3+(1-x)**-3
                refs=[ss0*(bb+16)/dd,3*ss0*(bb+16)*(x**-4-(1-x)**-4)/dd**2]
                for name, actual, ref in zip(["sigma","prime"],[surface_sigma(y,s,b),surface_sigma_prime(y,s,b)],refs):
                    rounded=float(ref)
                    err=abs(actual-rounded)
                    rel=err/abs(rounded) if abs(rounded)>=2.2250738585072014e-308 else None
                    ulps=err/math.ulp(rounded)
                    ok=math.isfinite(actual) and (rel<=5e-14 if rel is not None else ulps<=2)
                    rows.append(dict(function=name,y=y,b=b,sigma0=s,actual=actual,expected=rounded,relative_error=rel,ulp_error=ulps,ok=ok))
    for f in fs:
        ff=D.from_float(f)
        root=ctx.power(ff,D(1)/3) if f else D(0)
        ref=((2-3*root+ff)/5)
        rounded=float(ref)
        actual=spherical_phi(f)
        rel=abs(actual-rounded)/abs(rounded) if rounded else 0.
        rows.append(dict(function="phi",f=f,actual=actual,expected=rounded,relative_error=rel,ok=rel<=5e-14))
p=pathlib.Path(__file__).parent
src=pathlib.Path("python_frontend/compose_eos/scurto_interface.py")
result=dict(commit=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
 source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),python=platform.python_version(),
 decimal_precision=400,case_count=len(rows),failed=sum(not r["ok"] for r in rows),
 max_normal_relative_error=max(r["relative_error"] or 0 for r in rows),
 max_subnormal_ulp_error=max([r.get("ulp_error",0) for r in rows if r.get("relative_error") is None],default=0),
 rows=rows)
(p/"result.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({k:v for k,v in result.items() if k!="rows"},indent=2))
if result["failed"]:
    print(json.dumps([r for r in rows if not r["ok"]][:10],indent=2))
raise SystemExit(bool(result["failed"]))
