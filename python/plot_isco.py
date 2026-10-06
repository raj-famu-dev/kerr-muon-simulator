#!/usr/bin/env python3
"""Survival map for circular orbits: which radii are stable?

Reads the scan written by monograph_isco and draws, for prograde and
retrograde orbits, which radii survive (green) and which plunge (red),
with Bardeen's analytic ISCO radii marked for comparison.

Usage:
    python plot_isco.py isco_scan.csv --spin 0.9 --out isco_survival.png
"""
import argparse
import sys

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BG, FG, GRID = "#0b0720", "#e8e0ff", "#2a2250"
GREEN, RED, ORANGE, PURPLE = "#3ddc97", "#ff4d6d", "#ff9d2e", "#c04bff"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
        "axes.edgecolor": GRID, "axes.labelcolor": FG, "xtick.color": FG,
        "ytick.color": FG, "text.color": FG, "grid.color": GRID, "font.size": 11,
    })


def isco_analytic(a, prograde=True):
    """Bardeen's closed-form ISCO radius (units of M)."""
    z1 = 1 + (1 - a * a) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
    z2 = np.sqrt(3 * a * a + z1 * z1)
    root = np.sqrt((3 - z1) * (3 + z1 + 2 * z2))
    return 3 + z2 - root if prograde else 3 + z2 + root


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="isco_scan.csv from monograph_isco")
    p.add_argument("--spin", type=float, default=0.9, help="spin a used for the scan (default 0.9)")
    p.add_argument("--out", default=None, help="output image (default: next to the CSV)")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    args = p.parse_args()

    apply_style()
    df = pd.read_csv(args.csv)
    r = df["r"].to_numpy()
    width = float(np.median(np.diff(r)))
    a = args.spin
    r_plus = 1 + np.sqrt(1 - a * a)
    pro, ret = isco_analytic(a, True), isco_analytic(a, False)

    fig, ax = plt.subplots(figsize=(12, 4.2))
    rows = [("prograde_stable", 1.0, "Prograde\n(with the spin)"),
            ("retrograde_stable", 0.0, "Retrograde\n(against the spin)")]
    for col, y, _ in rows:
        colors = [GREEN if s else RED for s in df[col]]
        ax.bar(r, 0.7, width=width * 1.02, bottom=y - 0.35, color=colors, alpha=0.9, linewidth=0, antialiased=False)

    ax.axvspan(0, r_plus, color="black", zorder=3)
    ax.axvline(r_plus, color=ORANGE, lw=2, zorder=4)
    ax.text(r_plus - 0.05, 1.62, f"horizon\n{r_plus:.2f} M", ha="right", va="top", color=ORANGE, fontsize=10)

    for val, y, name in [(pro, 1.0, "prograde"), (ret, 0.0, "retrograde")]:
        ax.plot([val, val], [y - 0.45, y + 0.45], color="white", lw=2, zorder=5)
        ax.text(val, y + 0.5, f"ISCO {val:.2f} M", ha="center", va="bottom", fontsize=10, color="white")

    ax.set_yticks([1, 0])
    ax.set_yticklabels([lbl for _, _, lbl in rows])
    ax.set_xlim(r.min() - 0.1, r.max() + 0.1)
    ax.set_ylim(-0.6, 1.75)
    ax.set_xlabel("orbital radius r (units of M)")
    ax.set_title(f"Which circular orbits survive?  Kerr black hole, a = {a:g}", pad=12)
    ax.text(0.99, 0.04, "green = survives    red = plunges    white tick = Bardeen's analytic ISCO",
            transform=ax.transAxes, ha="right", fontsize=9, color=FG)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    out = args.out or str(args.csv).rsplit(".", 1)[0] + "_survival.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print(f"Prograde ISCO  {pro:.3f} M   Retrograde ISCO {ret:.3f} M   (difference {ret - pro:.2f} M)")
    print(f"Saved {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
