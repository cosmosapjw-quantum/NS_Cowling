from pathlib import Path
import subprocess,json,time,os
r=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
o=r/"data/low_density_crust_20261117/stellar/ddmex_comparison"
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python",str(r/"data/low_density_crust_20261117/stellar/finalize_mode_receipts.py")]
t=time.time()
with (o/"FINALIZE.stdout.log").open("w") as a,(o/"FINALIZE.stderr.log").open("w") as b:p=subprocess.run(argv,cwd=r,stdout=a,stderr=b)
(o/"FINALIZE_PROCESS.json").write_text(json.dumps(dict(argv=argv,actual_exit_code=p.returncode,seconds=time.time()-t),indent=2)+"\n")
print((o/"FINALIZE.stdout.log").read_text())
print((o/"FINALIZE.stderr.log").read_text())
print("ACTUAL_EXIT",p.returncode)
raise SystemExit(p.returncode)
