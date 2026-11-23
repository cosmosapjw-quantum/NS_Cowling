#!/usr/bin/env python3
"""Read the published matched CompOSE tables and export n, energy and pressure in physical units."""
import argparse
import csv
import io
import json
from pathlib import Path
import zipfile

from audit_remaining_reconstruction_matching import stored_text


def export(release, dest):
    count = 0
    with zipfile.ZipFile(dest, "w") as target:
        for directory in sorted(release.iterdir()):
            archives = directory/"PUBLIC_MATCHED_EOS.zip"
            if not archives.is_file():
                continue
            with zipfile.ZipFile(archives) as source, (directory/"ALL_CELLS.csv").open() as f:
                for row in csv.DictReader(f):
                    folder = row["candidate_directory"]
                    if not folder:
                        continue
                    if row["private"] == "True" or row["qualified"] != "True":
                        raise ValueError("Unqualified or private table in public archive")
                    def read(member):
                        return source.read(folder+"/"+member).decode().splitlines()
                    nb = read("eos.nb")
                    i0, i1 = int(nb[0]), int(nb[1])
                    density = [float(v) for v in nb[2:]]
                    therm = read("eos.thermo")
                    mn = float(therm[0].split()[0])
                    phases = {int(cols[1]):int(cols[3]) for line in read("eos.compo") if (cols:=line.split())}
                    if i1-i0+1 != len(density) or len(therm)-1 != len(density) or len(phases) != len(density):
                        raise ValueError("Incomplete CompOSE data")
                    receipt = json.loads(source.read(folder+"/MATCH.json"))["native"]
                    left = receipt.get("n1",receipt.get("n_low"))
                    right = receipt.get("n2",receipt.get("n_high"))
                    if not left or not right:
                        raise ValueError("Missing native transition bounds")
                    output = io.StringIO()
                    w = csv.writer(output,lineterminator="\n")
                    w.writerow(["nb_fm3","epsilon_MeV_fm3","pressure_MeV_fm3","Iphase","region"])
                    for idx, (n, line) in enumerate(zip(density,therm[1:]),i0):
                        c = line.split()
                        if int(c[1]) != idx:
                            raise ValueError("Density coordinate differs from thermodynamic serialization")
                        region = "crust" if n <= left else ("core" if n >= right else "interface")
                        w.writerow([repr(n),repr(n*(1.+float(c[9]))*mn),repr(n*float(c[3])),phases[idx],region])
                    stored_text(target,f"{directory.name}/{folder}.csv",output.getvalue())
                    count += 1
        stored_text(target,"README.txt","One file per public candidate that passed the retained crust audit.\n"
            "Join to ALL_CELLS.csv on reconstruction (directory) and census_cell_id (filename).\n"
            "The original matched EOS and native endpoint receipts are in PUBLIC_MATCHED_EOS.zip per model.\n"
            "At zero temperature epsilon = n * mn * (1+Q7), P=n*Q1. Region includes the respective source endpoints.\n"
            "Projected NV folders contain diagnostics subject to their source restrictions.\n")
    return count


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--release",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    print(export(a.release,a.out))
