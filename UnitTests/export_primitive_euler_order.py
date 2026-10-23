#!/usr/bin/env python3
"""Export spatial-order diagnostics from the focused Catch2 log; optionally plot them."""
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "Docs")
    parser.add_argument("--plot", action="store_true", help="Requires matplotlib; not needed to run the tests")
    args = parser.parse_args()
    lines = args.log.read_text().splitlines()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for prefix, name in (("ORDER_AUDIT,", "regions"), ("FACE_MOMENT,", "moments")):
        records = [line.removeprefix(prefix) for line in lines if line.startswith(prefix)]
        if len(records) < 2:
            raise ValueError("No complete diagnostics for " + prefix)
        (args.output_dir / f"primitive_euler_phase4_order_{name}.csv").write_text("\n".join(records) + "\n")
    if not args.plot:
        return
    import matplotlib.pyplot as plt

    with (args.output_dir / "primitive_euler_phase4_order_regions.csv").open() as stream:
        rows = list(csv.DictReader(stream))

    def curve(kind, variant, group):
        selected = [r for r in rows if r["kind"] == kind and r["variant"] == variant and r["group"] == group
                    and r["family"] == "regular" and r["gradient"] == "WEIGHTED_LEAST_SQUARES"]
        selected.sort(key=lambda r: int(r["N"]))
        return [int(r["N"]) for r in selected], [float(r["L2_group"]) for r in selected]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout="constrained")
    for variant, label in (("NODE_POINT:native", "Affine / NODE_POINT"),
                           ("MEAN_COEFF:native", "Affine / MEAN_COEFF")):
        n, error = curve("residual", variant, "all")
        axes[0].loglog(n, error, "o-", label=label)
    n, error = curve("residual", "NODE_POINT:native", "all")
    axes[0].loglog(n, [error[0] * (n[0] / x)**1.5 for x in n], "k:", label="Riferimento N^-1.5")
    axes[0].loglog(n, [error[0] * (n[0] / x)**2 for x in n], "k--", label="Riferimento N^-2")
    axes[0].set_title("L2 globale: operatore affine")
    for group, label in (("boundary", "Bordo fisico"), ("ring1", "Primo anello"), ("deep3", "Interno da tre strati")):
        n, error = curve("residual", "NODE_POINT:native", group)
        axes[1].loglog(n, error, "o-", label=label)
    axes[1].set_title("Errori regionali: affine / NODE_POINT")
    for ax in axes:
        ax.set_xlabel("N (suddivisioni per lato)")
        ax.set_ylabel("L2 pesata sui volumi duali")
        ax.set_xticks(n, [str(x) for x in n])
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8)
    target = args.output_dir / "primitive_euler_phase4_order.svg"
    fig.savefig(target)
    target.write_text("\n".join(line.rstrip() for line in target.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    main()
