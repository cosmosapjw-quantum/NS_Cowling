r"""Command-line interface for the ``compose_eos`` package (v3.1.0).

Supports:
  - ``compose-eos check <dir>`` (or ``compose-eos audit <dir>``)
  - ``compose-eos reconstruct <dir> -o <out>``
  - ``compose-eos rebuild <dir> -o <out>``
  - ``compose-eos match --crust <crust_dir> --core <core_dir> --out <out_dir> [--method all|...] [--tov]``
"""
from __future__ import annotations

import argparse
import os
import sys

from .eos_data import EOSData
from .io import read_compose, write_compose
from .consistency import check_consistency
from .reconstruct import fill_density_gap
from .matching import METHODS, match_crust_core, solve_star_at_target_mass


def _cmd_check(args) -> int:
    data = read_compose(args.input)
    print(data.summary())
    print()
    report = check_consistency(data, verbose=True)
    return 0 if report["overall"]["pass_"] else 2


def _cmd_reconstruct(args) -> int:
    data = read_compose(args.input)
    print(data.summary())
    print()
    print(f"filling density gap with {args.n_fill} points ...")
    rec = fill_density_gap(data, n_fill=args.n_fill)
    print(rec.summary())
    print()
    print("consistency of the reconstructed table:")
    report = check_consistency(rec, verbose=True)
    if not report["overall"]["pass_"]:
        print("\nrefusing to write a reconstructed EOS that fails hard consistency gates", file=sys.stderr)
        return 2
    write_compose(args.output, rec, copy_from=args.input)
    print(f"\nreconstructed table written to {args.output}")

    if args.figures:
        from .plotting import make_figures
        paths = make_figures(rec, args.figures, mr_path=args.mr)
        print(f"wrote {len(paths)} figures to {args.figures}")
    return 0


def _cmd_rebuild(args) -> int:
    from .reconstruct import reconstruct_master
    data = read_compose(args.input)
    print(data.summary())
    print()
    print("rebuilding from free energy as master function ...")
    rec = reconstruct_master(data, smooth=not args.no_smooth,
                             smooth_sigma=args.smooth_sigma)
    print(rec.summary())
    print()
    print("consistency of the rebuilt table:")
    report = check_consistency(rec, verbose=True)
    if not report["overall"]["pass_"]:
        print("\nrefusing to write a rebuilt EOS that fails hard consistency gates", file=sys.stderr)
        return 2
    write_compose(args.output, rec, copy_from=args.input)
    print(f"\nrebuilt table written to {args.output}")

    if args.figures:
        from .plotting import make_figures
        paths = make_figures(rec, args.figures, mr_path=args.mr)
        print(f"wrote {len(paths)} figures to {args.figures}")
    return 0


def _cmd_match(args) -> int:
    d_cr = read_compose(args.crust)
    d_co = read_compose(args.core)
    cr_name = args.crust_name or os.path.basename(os.path.normpath(args.crust)).replace("_fixed", "")
    co_name = args.core_name or os.path.basename(os.path.normpath(args.core)).replace("_fixed", "")
    methods = list(METHODS) if args.method == "all" else [args.method]

    for m in methods:
        eos_m, info = match_crust_core(
            d_cr, d_co, cr_name, co_name, method=m,
            n_bridge=args.n_bridge, core_compo_mode=args.core_compo_mode
        )
        rep = check_consistency(eos_m, verbose=False)
        if not rep["overall"]["pass_"] and not args.allow_failed:
            print(f"[{m}] hard consistency gate failed; refusing to write. Use --allow-failed only for an explicitly diagnostic export.", file=sys.stderr)
            return 2
        out_sub = os.path.join(args.out, f"{cr_name}__{co_name}__{m}") if len(methods) > 1 else args.out
        write_compose(out_sub, eos_m, copy_from=args.crust)
        msg = (
            f"[{m:10s}] construction={info['construction_kind']} | "
            f"dcs2_L={info['dcs2_left']:.2e}, dcs2_R={info['dcs2_right']:.2e} | "
            f"C2={rep['C2']['max_rel']:.2e}, overall={rep['overall']['status']}"
        )
        if args.tov:
            M14, R14, k2_14, Lam14, dR14, _ = solve_star_at_target_mass(
                eos_m, info["P1"], info["P2"], 80.0, 1.4
            )
            msg += f" | M={M14:.4f} Msun, R14={R14:.4f} km, k2={k2_14:.5f}, Lambda14={Lam14:.2f}, dR_cr={dR14:.4f} km"
        print(msg)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="compose-eos",
        description="Parse, check, reconstruct, match, and plot a CompOSE EOS table.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    for cmd_name in ("check", "audit"):
        pc = sub.add_parser(cmd_name, help="run thermodynamic consistency checks")
        pc.add_argument("input", help="CompOSE table directory")
        pc.set_defaults(func=_cmd_check)

    pr = sub.add_parser("reconstruct", help="fill the density gap and write a new table")
    pr.add_argument("input", help="source CompOSE table directory")
    pr.add_argument("-o", "--output", required=True, help="output directory")
    pr.add_argument("--n-fill", type=int, default=60, help="points to insert (default 60)")
    pr.add_argument("--figures", help="optional directory for diagnostic figures")
    pr.add_argument("--mr", help="optional path to eos.mr for the mass-radius figure")
    pr.set_defaults(func=_cmd_reconstruct)

    pb = sub.add_parser("rebuild",
                        help="rebuild from free energy as master function "
                             "(Gibbs-Duhem-exact, fixes composition & beta eq)")
    pb.add_argument("input", help="source CompOSE table directory")
    pb.add_argument("-o", "--output", required=True, help="output directory")
    pb.add_argument("--no-smooth", action="store_true",
                    help="disable F smoothing (keep shell wiggles)")
    pb.add_argument("--smooth-sigma", type=float, default=0.05,
                    help="F smoothing rms amplitude in MeV (default 0.05)")
    pb.add_argument("--figures", help="optional directory for diagnostic figures")
    pb.add_argument("--mr", help="optional path to eos.mr for the mass-radius figure")
    pb.set_defaults(func=_cmd_rebuild)

    pm = sub.add_parser("match", help="match a crust CompOSE table and a core CompOSE table")
    pm.add_argument("--crust", required=True, help="crust CompOSE directory")
    pm.add_argument("--core", required=True, help="core CompOSE directory")
    pm.add_argument("-o", "--out", required=True, help="output directory")
    pm.add_argument("--crust-name", default="", help="optional crust model label")
    pm.add_argument("--core-name", default="", help="optional core model label")
    pm.add_argument("--method", default="quintic_P", choices=list(METHODS) + ["all"],
                    help="matching scheme (default: quintic_P, or 'all')")
    pm.add_argument("--n-bridge", type=int, default=40, help="number of bridge points (default: 40)")
    pm.add_argument("--core-compo-mode", default="literature_beta_eq",
                    choices=["literature_beta_eq", "strict_thermo_only"],
                    help="core composition mode when eos.compo is absent on core input")
    pm.add_argument("--tov", action="store_true", help="compute 1.4 Msun TOV + tidal benchmarks")
    pm.set_defaults(func=_cmd_match)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
