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
    # E390-10B phase 1: full multidim tensor view (opt-in, unallocated by
    # default). q_full has shape (nT, nYq, 7, n), matching the Fortran
    # eos_data_t storage order q_full(nT,nYq,7,n) so select_slice can copy a
    # contiguous (7,n) slice per cell; q_extra_full has shape
    # (nadd_max, nT, nYq, n). An empty size-0 array means "not materialized";
    # cold-slice readers (the current EOS390 production path) leave these
    # empty so M3-M6 and the app commands keep working unchanged. select_slice
    # (api.select_slice) materializes a (iT,iYq) cell into the 1D Q*/q_extra
    # surface. nadd_max is the declared max extra-quantity count (0 for
    # Nadd=0 tables, 1 for the current cohort).
    nadd_max: int = 0
    q_full: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=float))
    q_extra_full: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=float))
    source_dir: str | None = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Initialize masks for directly constructed full-composition records.

        Programmatic EOSData constructors supply every composition array. If
        composition_available is true and no row masks were provided, those
        supplied values (including zeros) are explicit. Parser-created objects
        carry length-n masks and therefore retain their exact sparse coverage.
        """
        n = self.n
        presence_names = (
            "Yn_present", "Yp_present", "Ye_present", "Ymu_present",
            "Abar_present", "Zbar_present", "Yh_present",
        )
        state_names = presence_names + ("C3_applicable", "C4_applicable", "row_status")
        all_empty = all(np.asarray(getattr(self, name)).size == 0 for name in state_names)
        if not all_empty:
            return
        if self.composition_available:
            for name in presence_names:
                setattr(self, name, np.ones(n, dtype=bool))
            self.C3_applicable = self.Yn_present & self.Yp_present & self.Yh_present
            self.C4_applicable = (
                self.Yp_present & self.Abar_present & self.Zbar_present &
                self.Yh_present & self.Ye_present & self.Ymu_present & (self.Aav > 0.0)
            )
            self.row_status = np.where(
                self.C3_applicable & self.C4_applicable, 2, 1
            ).astype(np.int64)
        else:
            for name in presence_names:
                setattr(self, name, np.zeros(n, dtype=bool))
            self.C3_applicable = np.zeros(n, dtype=bool)
            self.C4_applicable = np.zeros(n, dtype=bool)
            self.row_status = np.full(
                n, 1 if self.composition_file_present else 0, dtype=np.int64
            )

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
        # E390-10B phase 3/4: when the CSV carries it_local and/or iy_local
        # columns, the rows span nT*nYq*n (T,Yq,density) triples. The 1D surface
        # (nb,Q*,q_extra,compo,masks) only materializes the selected (first)
        # T,Yq (it_local==1 and iy_local==1); the remaining T,Yq blocks feed the
        # q_full/q_extra_full tensors below.
        has_it_axis = "it_local" in header
        has_iyq_axis = "iy_local" in header
        if has_it_axis and has_iyq_axis:
            surf_rows = [r for r in rows if int(r["it_local"]) == 1 and int(r["iy_local"]) == 1]
        elif has_it_axis:
            surf_rows = [r for r in rows if int(r["it_local"]) == 1]
        elif has_iyq_axis:
            surf_rows = [r for r in rows if int(r["iy_local"]) == 1]
        else:
            surf_rows = rows
        # keep `rows` for tensor materialization; the col() helper below reads
        # from surf_rows (the 1D surface) so the existing field lookups stay
        # unchanged.
        orig_rows = rows
        rows = surf_rows
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
        eos = cls(
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
            # E390-10B phase 1/3: nadd_max + the axis lengths round-trip via the
            # namelist. q_full/q_extra_full tensors are materialized below from
            # the CSV (iy_local axis + q_extra_N columns).
            nadd_max=int(metadata.get("nadd_max", 0)),
            nT=int(metadata.get("naxis_t", 1)) or 1,
            nYq=int(metadata.get("naxis_yq", 1)) or 1,
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
        # E390-10B phase 2: collect the additional extra-quantity columns
        # (q_extra_2, q_extra_3, ...) the Fortran writer now emits when nadd_max
        # >= 2, and materialize them into the q_extra_full tensor
        # (nadd_max, nT, nYq, n). Slot 1 mirrors the 1D q_extra surface. Absent
        # columns default to 0 (backward compatible with Nadd in {0,1}).
        nadd_max_out = int(metadata.get("nadd_max", 0))
        nT = int(metadata.get("naxis_t", 1)) or 1
        nYq = int(metadata.get("naxis_yq", 1)) or 1
        import numpy as _np
        has_it_axis_qe = "it_local" in header
        has_iyq_axis_qe = "iy_local" in header
        n_density_qe = len(orig_rows) // (max(nT, 1) * max(nYq, 1)) if (has_it_axis_qe or has_iyq_axis_qe) else len(rows)
        if nadd_max_out >= 1:
            # E390-10B phase 3/4: when it_local/iy_local are present, the extra
            # columns span nT*nYq*n rows; collect per-(T,Yq) blocks into
            # q_extra_full. Otherwise (cold slice) the surface rows carry slot 1.
            qef = _np.zeros((nadd_max_out, nT, nYq, n_density_qe), dtype=float)
            if has_it_axis_qe and has_iyq_axis_qe:
                for it in range(1, nT + 1):
                    for iy in range(1, nYq + 1):
                        qef[0, it - 1, iy - 1, :] = _np.asarray(
                            [float(r.get("q_extra", 0.0) or 0.0)
                             for r in orig_rows
                             if int(r["it_local"]) == it and int(r["iy_local"]) == iy]
                        )
                for j in range(2, nadd_max_out + 1):
                    cname = f"q_extra_{j}"
                    if cname in header:
                        for it in range(1, nT + 1):
                            for iy in range(1, nYq + 1):
                                block = _np.asarray(
                                    [float(r[cname]) for r in orig_rows
                                     if int(r["it_local"]) == it and int(r["iy_local"]) == iy],
                                    dtype=float,
                                )
                                qef[j - 1, it - 1, iy - 1, :] = block
            elif has_iyq_axis_qe:
                # slot 1 (the q_extra column) for every Yq block
                for iy in range(1, nYq + 1):
                    qef[0, 0, iy - 1, :] = _np.asarray(
                        [float(r.get("q_extra", 0.0) or 0.0)
                         for r in orig_rows if int(r["iy_local"]) == iy]
                    )
                for j in range(2, nadd_max_out + 1):
                    cname = f"q_extra_{j}"
                    if cname in header:
                        for iy in range(1, nYq + 1):
                            block = _np.asarray(
                                [float(r[cname]) for r in orig_rows if int(r["iy_local"]) == iy],
                                dtype=float,
                            )
                            qef[j - 1, 0, iy - 1, :] = block
            elif has_it_axis_qe:
                for it in range(1, nT + 1):
                    qef[0, it - 1, 0, :] = _np.asarray(
                        [float(r.get("q_extra", 0.0) or 0.0)
                         for r in orig_rows if int(r["it_local"]) == it]
                    )
                for j in range(2, nadd_max_out + 1):
                    cname = f"q_extra_{j}"
                    if cname in header:
                        for it in range(1, nT + 1):
                            block = _np.asarray(
                                [float(r[cname]) for r in orig_rows if int(r["it_local"]) == it],
                                dtype=float,
                            )
                            qef[j - 1, it - 1, 0, :] = block
            else:
                qef[0, 0, 0, :] = _np.asarray([float(r.get("q_extra", 0.0) or 0.0) for r in rows])
                for j in range(2, nadd_max_out + 1):
                    cname = f"q_extra_{j}"
                    if cname in header:
                        qef[j - 1, 0, 0, :] = _np.asarray(
                            [float(r[cname]) for r in rows], dtype=float
                        )
            eos.q_extra_full = qef
        # E390-10B phase 3/4: materialize q_full(nT,nYq,7,n) when the CSV carries
        # an it_local and/or iy_local axis. The density count is len(rows)/(nT*nYq);
        # each (T,Yq) block holds n_density rows. Without it_local/iy_local the
        # tensor stays empty (legacy cold slice).
        has_it_axis = "it_local" in header
        has_iyq_axis = "iy_local" in header
        n_density = len(orig_rows) // (max(nT, 1) * max(nYq, 1)) if (has_it_axis or has_iyq_axis) else len(rows)
        if has_it_axis and has_iyq_axis and (nT > 1 or nYq > 1):
            qf = _np.zeros((nT, nYq, 7, n_density), dtype=float)
            for it in range(1, nT + 1):
                for iy in range(1, nYq + 1):
                    block = [r for r in orig_rows
                             if int(r["it_local"]) == it and int(r["iy_local"]) == iy]
                    for qi, qn in enumerate(("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7")):
                        qf[it - 1, iy - 1, qi, :] = _np.asarray(
                            [float(r[qn]) for r in block], dtype=float
                        )
            eos.q_full = qf
        elif has_iyq_axis and nYq > 1:
            qf = _np.zeros((nT, nYq, 7, n_density), dtype=float)
            for iy in range(1, nYq + 1):
                block = [r for r in orig_rows if int(r["iy_local"]) == iy]
                for qi, qn in enumerate(("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7")):
                    qf[0, iy - 1, qi, :] = _np.asarray(
                        [float(r[qn]) for r in block], dtype=float
                    )
            eos.q_full = qf
        elif has_it_axis and nT > 1:
            qf = _np.zeros((nT, nYq, 7, n_density), dtype=float)
            for it in range(1, nT + 1):
                block = [r for r in orig_rows if int(r["it_local"]) == it]
                for qi, qn in enumerate(("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7")):
                    qf[it - 1, 0, qi, :] = _np.asarray(
                        [float(r[qn]) for r in block], dtype=float
                    )
            eos.q_full = qf
        return eos

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
        # E390-10B phase 2: append one column per extra slot beyond the first.
        # Slot 1 is the existing q_extra; slots 2..nadd_max come from the
        # q_extra_full tensor (cold slice: indices nT-1, nYq-1). The tensor is
        # optional; absent slots are written as 0 and read back as 0.
        # E390-10B phase 3: when nYq>1, prepend an iy_local column and emit one
        # (Yq,density) row per Yq block (nYq*n rows total). Q1..Q7 for non-first
        # Yq blocks come from q_full; compo/mask columns are meaningful only for
        # the selected (first) Yq and zero-filled for the rest.
        nadd_max = int(self.nadd_max)
        nYq = int(self.nYq) or 1
        nT = int(self.nT) or 1
        qef = np.asarray(self.q_extra_full)
        qf = np.asarray(self.q_full)
        # q_full shape is (nT, nYq, 7, n); q_extra_full shape is (nadd_max, nT, nYq, n)
        has_iyq = nYq > 1 and qf.ndim == 4 and qf.shape[1] == nYq
        has_it = nT > 1 and qf.ndim == 4 and qf.shape[0] == nT
        extra_slot_cols = []
        for j in range(2, max(1, nadd_max) + 1):
            fields.append(f"q_extra_{j}")
            extra_slot_cols.append(j)  # resolved per-(T,Yq) below
        if has_it and has_iyq:
            fields = ["it_local", "iy_local"] + fields
        elif has_it:
            fields = ["it_local"] + fields
        elif has_iyq:
            fields = ["iy_local"] + fields
        arrays_spec = [self.nb, self.Q1, self.Q2, self.Q3, self.Q4, self.Q5, self.Q6, self.Q7,
                       self.q_extra, self.Iphase, self.Yn, self.Yp, self.Ye, self.Ymu,
                       self.Aav, self.Zav, self.Yh,
                       rs, yn_p, yp_p, ye_p, ymu_p, abar_p, zbar_p, yh_p, c3, c4]
        # integer columns (Iphase=idx 9->10 with it/iy_local prefix, row_status);
        # bool columns follow. Indices are computed relative to arrays_spec
        # (the it/iy_local prefix is written separately, not via arrays_spec).
        int_idx = {9, 17}
        bool_idx = {18, 19, 20, 21, 22, 23, 24, 25, 26}
        nrow_blocks_t = nT if has_it else 1
        nrow_blocks_y = nYq if has_iyq else 1
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            for it in range(1, nrow_blocks_t + 1):
                for iy in range(1, nrow_blocks_y + 1):
                    for k in range(n):
                        row = []
                        if has_it and has_iyq:
                            row.append(it); row.append(iy)
                        elif has_it:
                            row.append(it)
                        elif has_iyq:
                            row.append(iy)
                        is_selected = (it == 1 and iy == 1)
                        for j, a in enumerate(arrays_spec):
                            if j in int_idx:
                                row.append(int(a[k]) if is_selected else 0)
                            elif j in bool_idx:
                                row.append((1 if bool(a[k]) else 0) if is_selected else 0)
                            else:
                                # thermo Q1..Q7 (idx 1..7), nb (idx 0), q_extra (idx 8):
                                # for non-selected T,Yq blocks pull from q_full / q_extra_full.
                                if (has_it or has_iyq) and 1 <= j <= 7:
                                    row.append(format(float(qf[it - 1, iy - 1, j - 1, k]), ".17e"))
                                elif (has_it or has_iyq) and j == 8:
                                    if qef.ndim == 4 and qef.shape[0] >= 1:
                                        row.append(format(float(qef[0, it - 1, iy - 1, k]), ".17e"))
                                    else:
                                        row.append("0.00000000000000000e+00")
                                elif (has_it or has_iyq) and j == 0:
                                    row.append(format(float(a[k]), ".17e"))
                                else:
                                    row.append(format(float(a[k]), ".17e"))
                        # extra slots 2..nadd_max per-(T,Yq) from q_extra_full
                        for j in range(2, max(1, nadd_max) + 1):
                            if qef.ndim == 4 and qef.shape[0] >= j:
                                row.append(format(float(qef[j - 1, it - 1, iy - 1, k]), ".17e"))
                            else:
                                row.append("0.00000000000000000e+00")
                        writer.writerow(row)
        (directory / "eos_data.nml").write_text(
            "&metadata\n"
            f" mn={self.mn:.17e}\n mp={self.mp:.17e}\n il={self.Il}\n iT={self.iT}\n"
            f" iYq={self.iYq}\n inb_start={self.inb_start}\n nadd={self.q_extra_count}\n"
            f" naxis_t={int(self.nT)}\n naxis_yq={int(self.nYq)}\n"
            f" temperature_mev={self.temperature_mev:.17e}\n"
            f" charge_fraction_grid={self.charge_fraction_grid:.17e}\n"
            f" composition_file_present={'.true.' if self.composition_file_present else '.false.'}\n"
            f" composition_available={'.true.' if self.composition_available else '.false.'}\n"
            f" beta_equilibrated={'.true.' if self.beta_equilibrated else '.false.'}\n"
            f" nadd_max={int(self.nadd_max)}\n/\n",
            encoding="utf-8",
        )
        return directory

    def copy(self) -> "EOSData":
        names = ("nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7","q_extra","Iphase",
                 "Yn","Yp","Ye","Ymu","Aav","Zav","Yh",
                 "row_status","Yn_present","Yp_present","Ye_present","Ymu_present",
                 "Abar_present","Zbar_present","Yh_present","C3_applicable","C4_applicable",
                 "t_axis","yq_axis","q_full","q_extra_full")
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
            beta_equilibrated=self.beta_equilibrated, source_dir=self.source_dir,
            nT=self.nT, nYq=self.nYq, t_first=self.t_first, t_last=self.t_last,
            yq_first=self.yq_first, yq_last=self.yq_last, selected_iT=self.selected_iT,
            selected_iYq=self.selected_iYq, nadd_max=self.nadd_max)
