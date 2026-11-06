"""End-to-end example: reconstruct EOS 217 and verify consistency.

Run with::

    python reconstruct_eos217.py /path/to/eos217 /path/to/output

If no arguments are given it falls back to the bundled paths used during
development.
"""
import sys

from compose_eos import read_compose, check_consistency, fill_density_gap, write_compose, make_figures


def main(src="eos217", out="eos217_reconstructed", figs="figures"):
    data = read_compose(src)
    print(data.summary())
    print("\n--- original table ---")
    check_consistency(data)

    rec = fill_density_gap(data, n_fill=60)
    print("\n--- reconstructed table ---")
    check_consistency(rec)

    write_compose(out, rec, copy_from=src)
    print(f"\nreconstructed table written to {out}")
    make_figures(rec, figs, mr_path=f"{src}/eos.mr")
    print(f"figures written to {figs}")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(*args) if args else main()
