"""Python-side EOSData container; all EOS transformations are delegated to Fortran."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import csv
import re
from typing import Any
import numpy as np

def _read_metadata(path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, raw = line.strip().rstrip(",").split("=", 1)
        key = key.strip().lower()
        value = raw.strip()
        logical_value = value.lower()
        if logical_value in (".true.", "true", "t"):
            result[key] = True
        elif logical_value in (".false.", "false", "f"):
            result[key] = False
        elif re.fullmatch(r"[+-]?\d+", value):
            result[key] = int(value)
        else:
            try:
                result[key] = float(value.replace("D", "E").replace("d", "e"))
            except ValueError:
                result[key] = value.strip("'\"")
    return result

@dataclass
class EOSData:
    nb: np.ndarray
    mn: float
    mp: float
    Q1: np.ndarray
    Q2: np.ndarray
    Q3: np.ndarray
    Q4: np.ndarray
    Q5: np.ndarray
    Q6: np.ndarray
    Q7: np.ndarray
    q_extra: np.ndarray
    Iphase: np.ndarray
    Yn: np.ndarray
    Yp: np.ndarray
    Ye: np.ndarray
    Ymu: np.ndarray
    Aav: np.ndarray
    Zav: np.ndarray
    Yh: np.ndarray
    Il: int = 1
    iT: int = 1
    iYq: int = 1
    inb_start: int = 1
    temperature_mev: float = 0.0
    charge_fraction_grid: float = 0.0
    q_extra_count: int = 0
    composition_file_present: bool = False
    composition_available: bool = False
    beta_equilibrated: bool = False
    source_dir: str | None = field(default=None, repr=False, compare=False)

    @property
    def n(self) -> int:
        return int(self.nb.size)
    @property
    def p(self) -> np.ndarray:
        return self.Q1 * self.nb
    @property
    def mub(self) -> np.ndarray:
        return (self.Q3 + 1.0) * self.mn
    @property
    def muq(self) -> np.ndarray:
        return self.Q4 * self.mn
    @property
    def mul(self) -> np.ndarray:
        return self.Q5 * self.mn
    @property
    def mue(self) -> np.ndarray:
        return -self.muq
    @property
    def F(self) -> np.ndarray:
        return (self.Q6 + 1.0) * self.mn
    @property
    def f(self) -> np.ndarray:
        return self.nb * self.F
    @property
    def eps(self) -> np.ndarray:
        return self.nb * (self.Q7 + 1.0) * self.mn
    @property
    def entropy(self) -> np.ndarray:
        return self.Q2

    @classmethod
    def from_internal_directory(cls, directory: str | Path, source_dir: str | None = None) -> "EOSData":
        directory = Path(directory)
        metadata = _read_metadata(directory / "eos_data.nml")
        with (directory / "eos_data.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            raise ValueError("Fortran backend produced an empty EOSData CSV")
        def col(name: str, dtype=float):
            return np.asarray([dtype(row[name]) for row in rows])
        return cls(
            nb=col("nb"), mn=float(metadata["mn"]), mp=float(metadata["mp"]),
            Q1=col("Q1"), Q2=col("Q2"), Q3=col("Q3"), Q4=col("Q4"),
            Q5=col("Q5"), Q6=col("Q6"), Q7=col("Q7"), q_extra=col("q_extra"),
            Iphase=col("Iphase", int), Yn=col("Yn"), Yp=col("Yp"), Ye=col("Ye"),
            Ymu=col("Ymu"), Aav=col("Abar"), Zav=col("Zbar"), Yh=col("Yh"),
            Il=int(metadata["il"]), iT=int(metadata["it"]), iYq=int(metadata["iyq"]),
            inb_start=int(metadata["inb_start"]), temperature_mev=float(metadata["temperature_mev"]),
            charge_fraction_grid=float(metadata["charge_fraction_grid"]), q_extra_count=int(metadata["nadd"]),
            composition_file_present=bool(metadata["composition_file_present"]),
            composition_available=bool(metadata["composition_available"]),
            beta_equilibrated=bool(metadata["beta_equilibrated"]), source_dir=source_dir,
        )

    def write_internal_directory(self, directory: str | Path) -> Path:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        csv_path = directory / "eos_data.csv"
        fields = ["nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra",
                  "Iphase", "Yn", "Yp", "Ye", "Ymu", "Abar", "Zbar", "Yh"]
        arrays = [self.nb, self.Q1, self.Q2, self.Q3, self.Q4, self.Q5, self.Q6, self.Q7,
                  self.q_extra, self.Iphase, self.Yn, self.Yp, self.Ye, self.Ymu,
                  self.Aav, self.Zav, self.Yh]
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            for i in range(self.n):
                writer.writerow([format(float(a[i]), ".17e") if j != 9 else int(a[i])
                                 for j, a in enumerate(arrays)])
        (directory / "eos_data.nml").write_text(
            "&metadata\n"
            f" mn={self.mn:.17e}\n mp={self.mp:.17e}\n il={self.Il}\n iT={self.iT}\n"
            f" iYq={self.iYq}\n inb_start={self.inb_start}\n nadd={self.q_extra_count}\n"
            f" temperature_mev={self.temperature_mev:.17e}\n"
            f" charge_fraction_grid={self.charge_fraction_grid:.17e}\n"
            f" composition_file_present={'.true.' if self.composition_file_present else '.false.'}\n"
            f" composition_available={'.true.' if self.composition_available else '.false.'}\n"
            f" beta_equilibrated={'.true.' if self.beta_equilibrated else '.false.'}\n/\n",
            encoding="utf-8",
        )
        return directory

    def copy(self) -> "EOSData":
        values = {name: np.array(getattr(self, name), copy=True) for name in
                  ("nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7","q_extra","Iphase",
                   "Yn","Yp","Ye","Ymu","Aav","Zav","Yh")}
        return EOSData(**values, mn=self.mn, mp=self.mp, Il=self.Il, iT=self.iT, iYq=self.iYq,
            inb_start=self.inb_start, temperature_mev=self.temperature_mev,
            charge_fraction_grid=self.charge_fraction_grid, q_extra_count=self.q_extra_count,
            composition_file_present=self.composition_file_present,
            composition_available=self.composition_available,
            beta_equilibrated=self.beta_equilibrated, source_dir=self.source_dir)
