#!/usr/bin/env python3
"""Regrid existing native N=128 TF branch seeds to N=256 input for the native driver.

Python only changes radial mesh sampling. EOS physics and minimization remain
in blocker_tf_common_tangent_polish.f90.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("FSU2R", "TM1e"))
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    raw = args.source.read_text().splitlines()
    n0 = 128
    n1 = 256
    cursor = 0
    output: list[str] = []
    source_x = (np.arange(n0) + 0.5) / n0
    target_x = (np.arange(n1) + 0.5) / n1
    for _ in range(2):
        header = raw[cursor].split()
        cursor += 1
        if int(header[0]) not in (1, 2, 3):
            raise ValueError(f"expected 1D branch seed, got {header[0]}")
        output.append(" ".join(header))
        fields = raw[cursor].split()
        cursor += 1
        if len(fields) != n0:
            raise ValueError(f"expected {n0} values, got {len(fields)}")
        values = np.asarray(fields, dtype=float)
        output.extend(f"{x:.16e}" for x in np.interp(target_x, source_x, values))
        fields = raw[cursor].split()
        cursor += 1
        if len(fields) != n0:
            raise ValueError(f"expected {n0} values, got {len(fields)}")
        values = np.asarray(fields, dtype=float)
        output.extend(f"{x:.16e}" for x in np.interp(target_x, source_x, values))
    if cursor != len(raw):
        raise ValueError(f"unexpected trailing input lines: {len(raw)-cursor}")
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    args.destination.write_text("\n".join(output) + "\n")
    print(f"model={args.model} n_source={n0} n_target={n1}")
    print(f"source_sha256={sha(args.source)}")
    print(f"destination_sha256={sha(args.destination)}")


if __name__ == "__main__":
    main()
