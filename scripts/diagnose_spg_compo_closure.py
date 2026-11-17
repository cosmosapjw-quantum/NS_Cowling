"""Diagnose the SPG_M1 eos.compo against the CompOSE C3/C4 closure identities.

CompOSE standard:
  C3 (baryon closure):  Yn + Yp + Yh = 1   where Yh = A*Yi (heavy-cluster baryon fraction)
  C4 (charge closure):  Yp + (Z/A)*Yh = Ye + Ymu

The provider eos.compo row tuple is (iN, iN_id, Yi_gas, ...) + (A, Z, Yi_cluster).
The doc says: "Do not recycle the non-conformant provider eos.compo group tuple
as a physics target" -- meaning the raw (A,Z,Yi) does not by itself satisfy the
CompOSE group-abundance standard. This script measures the gap.
"""
import numpy as np
from pathlib import Path

RAW = Path("inputs/raw_compose")


def parse_compo(model):
    rows = []
    with open(RAW / model / "eos.compo") as f:
        for line in f:
            t = line.split()
            if len(t) < 18:
                continue
            it, inb, iyq = int(t[0]), int(t[1]), int(t[2])
            iphase = int(t[4])
            m = int(t[5])  # length of the M-block that follows
            # the M-block layout (CompOSE compo v1):
            #   for each particle/group: (id, Yi, present)
            # neutron(11), proton(12), electron(13), muon(14), cluster(1)
            # then trailing (A, Z, Yi_cluster)
            # The actual column layout from the SPG rows:
            # col6=M=10, then 10 values: [id_n(11), Yi_n, pres_n? , id_p? ...]
            # We observed: col6=10, col7..15 = mixed id/Yi/present, col16=A, col17=Z, col18=Yi_cluster
            # Parse robustly: neutron block (id=11), proton block (id=12?), cluster block
            # From row: col5=4(Iphase), col6=10(M=10 numbers follow)
            block = t[6:6 + m]
            # cluster is the last triple in the block: (A, Z, Yi) -- but block is 10 wide
            # Actually col16=A, col17=Z, col18=Yi are AFTER the M-block? No: M=10 means
            # 10 values from col7..col16. But A is at col16. Let's just read by known position.
            # The observed stable layout (all SPG_M1 rows): 18 fields.
            # col1..4 = key+cluster_idx, col5=Iphase, col6=M=10,
            # col7,8,9  = (id=11, Yi_neutron, present)
            # col10,11,12 = (id=?, Yi_proton, present)
            # col13,14,15 = (id=1, Yi_cluster_marker, present)
            # col16=A, col17=Z, col18=Yi_cluster
            yi_n = float(t[8])
            yi_p = float(t[11])
            yi_cl = float(t[18 - 1])  # col18
            a_cl = float(t[16])
            z_cl = float(t[17])
            rows.append(dict(inb=inb, iphase=iphase, yi_n=yi_n, yi_p=yi_p,
                             a=a_cl, z=z_cl, yi_cl=yi_cl))
    return rows


def closure_residuals(rows, ye_data=None):
    """Compute C3 and C4 residuals. ye_data[i] = (Ye, Ymu) per row if available."""
    c3 = []
    c4 = []
    for i, r in enumerate(rows):
        yh = r["a"] * r["yi_cl"]  # heavy-cluster baryon fraction (Yh = A*Yi)
        c3_resid = r["yi_n"] + r["yi_p"] + yh - 1.0
        c3.append(c3_resid)
        if ye_data is not None and i < len(ye_data):
            ye, ymu = ye_data[i]
            # charge closure: Yp + (Z/A)*Yh = Ye + Ymu
            c4_resid = r["yi_p"] + (r["z"] / r["a"]) * yh - (ye + ymu)
            c4.append(c4_resid)
    return np.array(c3), np.array(c4) if c4 else None


for model in ["SPG_M1"]:
    rows = parse_compo(model)
    print(f"=== {model}: {len(rows)} compo rows ===")
    print(f"  Iphase values: {sorted(set(r['iphase'] for r in rows))}")
    print(f"  A range: [{min(r['a'] for r in rows):.3f}, {max(r['a'] for r in rows):.3f}]")
    print(f"  Z range: [{min(r['z'] for r in rows):.3f}, {max(r['z'] for r in rows):.3f}]")
    print(f"  Yi_n range: [{min(r['yi_n'] for r in rows):.5f}, {max(r['yi_n'] for r in rows):.5f}]")
    print(f"  Yi_p range: [{min(r['yi_p'] for r in rows):.5f}, {max(r['yi_p'] for r in rows):.5f}]")
    print(f"  Yi_cl range: [{min(r['yi_cl'] for r in rows):.5f}, {max(r['yi_cl'] for r in rows):.5f}]")
    c3, _ = closure_residuals(rows)
    yh = np.array([r["a"] * r["yi_cl"] for r in rows])
    print(f"\n  C3 (Yn+Yp+Yh=1) residuals: max|r|={np.max(np.abs(c3)):.3e}, mean={np.mean(c3):.3e}")
    print(f"  Yh = A*Yi_cl range: [{yh.min():.4f}, {yh.max():.4f}]")
    print(f"  Yn+Yp+Yh range: [{(np.array([r['yi_n'] for r in rows])+np.array([r['yi_p'] for r in rows])+yh).min():.6f}, {(np.array([r['yi_n'] for r in rows])+np.array([r['yi_p'] for r in rows])+yh).max():.6f}]")
    # show first/middle/last row detail
    for idx in [0, len(rows)//2, len(rows)-1]:
        r = rows[idx]
        yh_r = r["a"] * r["yi_cl"]
        print(f"  row {idx}: A={r['a']:.3f} Z={r['z']:.3f} Yi_n={r['yi_n']:.5f} Yi_p={r['yi_p']:.5f} Yi_cl={r['yi_cl']:.5f} Yh={yh_r:.5f} C3={r['yi_n']+r['yi_p']+yh_r-1:.3e}")
