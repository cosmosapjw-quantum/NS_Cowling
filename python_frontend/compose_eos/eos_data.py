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
    # M2A per-row presence masks: row_status in {0=NOT_AVAILABLE, 1=PARTIAL,
    # 2=CHECKED, 3=FAIL}; *_present is per-row logical (a written zero is True,
    # a missing entry is False); c3/c4_applicable per-row logical.
    row_status: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    Yn_present: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    Yp_present: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    Ye_present: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    Ymu_present: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    Abar_present: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    Zbar_present: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    Yh_present: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    C3_applicable: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    C4_applicable: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    # M10 generic CompOSE axes metadata (scalar; axis values stored separately).
    nT: int = 0
    nYq: int = 0
    t_first: int = 0
    t_last: int = -1
    yq_first: int = 0
    yq_last: int = -1
    selected_iT: int = 0
    selected_iYq: int = 0
    t_axis: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=float))
    yq_axis: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=float))
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
            reader = csv.DictReader(stream)
            rows = list(reader)
            header = reader.fieldnames or []
        if not rows:
            raise ValueError("Fortran backend produced an empty EOSData CSV")
        def col(name: str, dtype=float, default=None):
            if name not in header:
                if default is None:
                    raise ValueError(f"missing required column {name}")
                n = len(rows)
                if dtype is int:
                    return np.asarray([default] * n, dtype=np.int64)
                if dtype is bool:
                    return np.asarray([default] * n, dtype=bool)
                return np.asarray([default] * n, dtype=dtype)
            if dtype is bool:
                # the fortran writes 0/1 integers; convert via int to avoid
                # bool('0') being True (non-empty string).
                return np.asarray([bool(int(row[name])) for row in rows], dtype=bool)
            return np.asarray([dtype(row[name]) for row in rows])
        has_masks = "row_status" in header
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
            beta_equilibrated=bool(metadata["beta_equilibrated"]),
            # M2A masks: default to 0 / False if the backend wrote a legacy 17-col CSV.
            # NOTE: the Fortran CSV header uses lowercase column names (row_status,
            # yn_present, ...) which differ from the dataclass field names. The col
            # helper looks up by exact name, so use the lowercase header names here.
            row_status=col("row_status", int, default=0),
            Yn_present=col("yn_present", bool, default=False),
            Yp_present=col("yp_present", bool, default=False),
            Ye_present=col("ye_present", bool, default=False),
            Ymu_present=col("ymu_present", bool, default=False),
            Abar_present=col("abar_present", bool, default=False),
            Zbar_present=col("zbar_present", bool, default=False),
            Yh_present=col("yh_present", bool, default=False),
            C3_applicable=col("c3_applicable", bool, default=False),
            C4_applicable=col("c4_applicable", bool, default=False),
            source_dir=source_dir,
        )

    def write_internal_directory(self, directory: str | Path) -> Path:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        csv_path = directory / "eos_data.csv"
        n = self.n
        # M2A masks: if a caller constructed an EOSData directly (e.g. reconstruct_master)
        # the mask arrays may be empty. Derive safe defaults from the existing
        # composition arrays so legacy-constructed objects still round-trip:
        #   row_status: 2 (CHECKED) if composition_available else 0 (NOT_AVAILABLE)
        #   *_present: True where the corresponding array has a nonzero/positive entry
        #   C3/C4: not re-derived here (set by the parser); default False.
        def _mask_arr(arr, default, name):
            if arr is None:
                arr = np.zeros(0, dtype=bool)
            a = np.asarray(arr)
            if a.size == 0:
                # No composition sidecar means all rows are NOT_AVAILABLE. Once a
                # sidecar exists, missing row masks are ambiguous: they cannot be
                # reconstructed from numerical zeros, so fail closed.
                if self.composition_file_present and name.endswith("_present"):
                    raise ValueError(f"missing row-level presence mask {name}; refusing lossy export")
                return np.full(n, default)
            if a.shape != (n,):
                raise ValueError(f"row-level mask {name} has shape {a.shape}, expected ({n},)")
            return a

        rs = _mask_arr(getattr(self, "row_status", None), 0, "row_status").astype(int)
        def _bool_from(arr, default, name):
            return _mask_arr(arr, default, name).astype(bool)
        yn_p = _bool_from(getattr(self, "Yn_present", None), False, "Yn_present")
        yp_p = _bool_from(getattr(self, "Yp_present", None), False, "Yp_present")
        ye_p = _bool_from(getattr(self, "Ye_present", None), False, "Ye_present")
        ymu_p = _bool_from(getattr(self, "Ymu_present", None), False, "Ymu_present")
        abar_p = _bool_from(getattr(self, "Abar_present", None), False, "Abar_present")
        zbar_p = _bool_from(getattr(self, "Zbar_present", None), False, "Zbar_present")
        yh_p = _bool_from(getattr(self, "Yh_present", None), False, "Yh_present")
        c3 = _bool_from(getattr(self, "C3_applicable", None), False, "C3_applicable")
        c4 = _bool_from(getattr(self, "C4_applicable", None), False, "C4_applicable")
        fields = ["nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra",
                  "Iphase", "Yn", "Yp", "Ye", "Ymu", "Abar", "Zbar", "Yh",
                  # M2A per-row masks
                  "row_status", "Yn_present", "Yp_present", "Ye_present", "Ymu_present",
                  "Abar_present", "Zbar_present", "Yh_present",
                  "C3_applicable", "C4_applicable"]
        arrays = [self.nb, self.Q1, self.Q2, self.Q3, self.Q4, self.Q5, self.Q6, self.Q7,
                  self.q_extra, self.Iphase, self.Yn, self.Yp, self.Ye, self.Ymu,
                  self.Aav, self.Zav, self.Yh,
                  rs, yn_p, yp_p, ye_p, ymu_p, abar_p, zbar_p, yh_p, c3, c4]
        # integer columns (Iphase=10, row_status=18); bool columns are 19..26
        int_idx = {9, 17}
        bool_idx = {18, 19, 20, 21, 22, 23, 24, 25, 26}
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            for i in range(n):
                row = []
                for j, a in enumerate(arrays):
                    if j in int_idx:
                        row.append(int(a[i]))
                    elif j in bool_idx:
                        row.append(1 if bool(a[i]) else 0)
                    else:
                        row.append(format(float(a[i]), ".17e"))
                writer.writerow(row)
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
        names = ("nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7","q_extra","Iphase",
                 "Yn","Yp","Ye","Ymu","Aav","Zav","Yh",
                 "row_status","Yn_present","Yp_present","Ye_present","Ymu_present",
                 "Abar_present","Zbar_present","Yh_present","C3_applicable","C4_applicable")
        values = {}
        for name in names:
            a = getattr(self, name, None)
            if a is None:
                values[name] = np.zeros(0, dtype=np.int64 if name == "row_status" else bool)
            else:
                values[name] = np.array(a, copy=True)
        return EOSData(**values, mn=self.mn, mp=self.mp, Il=self.Il, iT=self.iT, iYq=self.iYq,
            inb_start=self.inb_start, temperature_mev=self.temperature_mev,
            charge_fraction_grid=self.charge_fraction_grid, q_extra_count=self.q_extra_count,
            composition_file_present=self.composition_file_present,
            composition_available=self.composition_available,
            beta_equilibrated=self.beta_equilibrated, source_dir=self.source_dir)
