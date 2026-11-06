"""Data model for a one-dimensional CompOSE equation-of-state table.

The :class:`EOSData` container holds the raw CompOSE quantities (the ``Q``
array and the composition) together with the physical thermodynamic variables
derived from them.  At :math:`T=0` and fixed charge fraction the table is
fully described by the baryon density grid ``nb``; every other quantity is a
function of ``nb``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class EOSData:
    r"""A one-dimensional, :math:`T=0` CompOSE equation-of-state table.

    Parameters
    ----------
    nb : ndarray
        Baryon-number density grid in fm\ :sup:`-3`.
    mn, mp : float
        Nucleon masses used by the table (MeV).  For DD2 ``mn = mp = 939.0``.
    Q1..Q7 : ndarray
        The seven CompOSE thermodynamic quantities (dimensionless, see the
        module docstring of :mod:`compose_eos.io` for their definitions).
    q_extra : ndarray
        Additional per-point quantities (here the pasta-geometry flag).
    Iphase : ndarray, int
        Phase index per point (0 outer crust, 7 inner crust, 6 core).
    Yn, Yp, Ye, Ymu : ndarray
        Particle fractions of neutrons, protons, electrons and muons.
    Aav, Zav, Yh : ndarray
        Mean baryon number, charge number, and *internal baryonic fraction* of
        the represented nuclear group/cell.  CompOSE stores the group number
        abundance Y_I; this package uses Yh = Aav * Y_I internally.
    Il : int
        CompOSE lepton-inclusion flag from the first eos.thermo line.
    iT, iYq, temperature_mev, charge_fraction_grid :
        The selected one-dimensional CompOSE slice coordinates.
    q_extra_count : Optional[int]
        Exact CompOSE Nadd; None means infer from q_extra when constructing
        new data by hand.  Existing input files set this explicitly.
    composition_available : bool
        Whether an eos.compo table existed.  Absence is valid in CompOSE.
    beta_equilibrated : bool
        Whether the selected cold slice is assumed to be beta equilibrated.
    """

    nb: np.ndarray
    mn: float = 939.0
    mp: float = 939.0
    Q1: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Q2: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Q3: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Q4: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Q5: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Q6: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Q7: np.ndarray = field(default_factory=lambda: np.zeros(0))
    q_extra: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Iphase: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    Yn: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Yp: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Ye: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Ymu: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Aav: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Zav: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Yh: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # CompOSE/source metadata retained across read -> reconstruction -> write.
    # Yh is an INTERNAL baryon fraction (Abar * CompOSE Y_I), never the raw
    # CompOSE nuclear number abundance.  io.py performs the conversion.
    Il: int = 1
    iT: int = 0
    iYq: int = 1
    inb_start: int = 1
    temperature_mev: float = 0.0
    charge_fraction_grid: float = 0.0
    q_extra_count: Optional[int] = None
    composition_available: bool = True
    beta_equilibrated: bool = True

    # ---- derived physical quantities (computed lazily) ----------------------
    @property
    def n(self) -> int:
        """Number of grid points."""
        return len(self.nb)

    @property
    def p(self) -> np.ndarray:
        r"""Pressure :math:`p = Q_1 n_b` in MeV fm\ :sup:`-3`."""
        return self.Q1 * self.nb

    @property
    def mub(self) -> np.ndarray:
        r"""Baryon chemical potential :math:`\\mu_b = (Q_3+1) m_n` in MeV."""
        return (self.Q3 + 1.0) * self.mn

    @property
    def muq(self) -> np.ndarray:
        r"""Charge chemical potential :math:`\\mu_q = Q_4 m_n` in MeV."""
        return self.Q4 * self.mn

    @property
    def mul(self) -> np.ndarray:
        r"""Effective lepton chemical potential :math:`\\mu_l = Q_5 m_n` in MeV."""
        return self.Q5 * self.mn

    @property
    def mue(self) -> np.ndarray:
        r"""Electron chemical potential :math:`\\mu_e = -\\mu_q` (beta eq.) in MeV."""
        return -self.muq

    @property
    def F(self) -> np.ndarray:
        """Free energy per baryon :math:`F = (Q_6+1) m_n` in MeV."""
        return (self.Q6 + 1.0) * self.mn

    @property
    def f(self) -> np.ndarray:
        r"""Free-energy density :math:`f = n_b F` in MeV fm\ :sup:`-3`."""
        return self.nb * self.F

    @property
    def eps(self) -> np.ndarray:
        r"""Energy density :math:`\\varepsilon = n_b (Q_7+1) m_n` in MeV fm\ :sup:`-3`."""
        return self.nb * (self.Q7 + 1.0) * self.mn

    @property
    def entropy(self) -> np.ndarray:
        """Entropy per baryon :math:`s = Q_2` (zero at :math:`T=0`)."""
        return self.Q2

    @property
    def ZoverA(self) -> np.ndarray:
        """Cell charge-to-baryon ratio :math:`Z/A` (zero where no nucleus)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.Aav > 0, self.Zav / self.Aav, 0.0)

    @property
    def Yn_eff(self) -> np.ndarray:
        r"""Effective free-neutron fraction from baryon conservation.

        Computed as ``1 - Yp - Yh``.  For tables that store the free-neutron
        fraction explicitly (e.g.\ EOS 134) this reproduces the stored value
        where the table is consistent; for tables that store it implicitly as
        ``Yn = 0`` in the pasta phase (e.g.\ EOS 217) it reconstructs the
        neutron-gas fraction so that ``Yn + Yp + Yh = 1`` everywhere.
        """
        return 1.0 - self.Yp - self.Yh

    # ---- phase masks --------------------------------------------------------
    # Convention-agnostic masks: a point "has a nucleus" when A > 0, regardless
    # of the integer Iphase code (EOS 217 uses 0/7/6; EOS 134 uses 4/3).
    @property
    def has_nucleus(self) -> np.ndarray:
        """Points where a nuclear cluster / cell is present (``A > 0``)."""
        return self.Aav > 0

    @property
    def is_homogeneous(self) -> np.ndarray:
        """Points with no nucleus (homogeneous core, ``A == 0``)."""
        return self.Aav == 0

    @property
    def is_outer(self) -> np.ndarray:
        """Legacy convenience mask only; phase ID 0 is not universal."""
        return self.Iphase == 0

    @property
    def is_inner(self) -> np.ndarray:
        """Legacy convenience mask only; phase ID 7 is not universal."""
        return self.Iphase == 7

    @property
    def is_core(self) -> np.ndarray:
        """Legacy convenience mask only; phase ID 6 is not universal."""
        return self.Iphase == 6

    # ---- summary ------------------------------------------------------------
    def summary(self) -> str:
        """Return a short human-readable summary of the table."""
        lines = [
            f"EOSData: {self.n} points, mn = {self.mn:.3f} MeV",
            f"  nb range: {self.nb.min():.3e} .. {self.nb.max():.3e} fm^-3",
            f"  p  range: {self.p.min():.3e} .. {self.p.max():.3e} MeV/fm^3",
            f"  mu_b range: {self.mub.min():.3f} .. {self.mub.max():.3f} MeV",
            f"  max |mu_l| = {np.max(np.abs(self.mul)):.3e} MeV",
            f"  max c_s^2 = {self._cs2_max():.4f}",
            f"  CompOSE slice: T={self.temperature_mev:g} MeV, Yq-grid={self.charge_fraction_grid:g}, Il={self.Il}",
            f"  composition: {'present' if self.composition_available else 'absent'}; Nadd={self.q_extra_count}",
        ]
        # CompOSE Iphase identifiers are provider/model-specific; never attach
        # universal physical labels to integer codes.
        if self.composition_available and self.Iphase.size == self.n:
            for ip in np.unique(self.Iphase):
                m = self.Iphase == ip
                lines.append(f"  phase ID {int(ip)}: {m.sum()} pts, "
                             f"nb {self.nb[m].min():.3e}..{self.nb[m].max():.3e}")
        return "\n".join(lines)

    def _cs2_max(self) -> float:
        deps = np.gradient(self.eps, self.nb)
        dP = np.gradient(self.p, self.nb)
        cs2 = np.where(np.abs(deps) > 0, dP / deps, 0.0)
        return float(np.nanmax(cs2))
