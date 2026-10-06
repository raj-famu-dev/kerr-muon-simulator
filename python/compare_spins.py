#!/usr/bin/env python3
"""Compare plunges side by side: a = 0, prograde and retrograde.

The muon always lives 2.2 microseconds of its own time. What changes with the
black hole's spin is WHERE along its path that happens.

Usage:
    python compare_spins.py schwarzschild.csv prograde.csv retrograde.csv \\
        --labels "a = 0" "a = 0.9 prograde" "a = 0.9 retrograde" --out compare_spins.png
"""
import argparse
import os
import sys

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle

BG, FG, GRID = "#0b0720", "#e8e0ff", "#2a2250"
PALETTE = ["#aab4ff", "#6ee7ff", "#ff9d2e", "#3ddc97", "#ff4d6d", "#c04bff"]
ORANGE = "#ff9d2e"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
        "axes.edgecolor": GRID, "axes.labelcolor": FG, "xtick.color": FG,
        "ytick.color": FG, "text.color": FG, "grid.color": GRID, "font.size": 11,
    })


def summarize(df):
    dead = df[df["muon_alive"] == 0]
    d = dead.iloc[0] if len(dead) else None
    return {
        "tau_end_us": float(df["tau_us"].iloc[-1]),
        "t_end_us": float(df["t_us"].iloc[-1]),
        "turns": float(df["phi"].iloc[-1]) / (2 * np.pi),
        "r_plus": float(df["r"].iloc[-1]) if float(df["dt_dtau"].iloc[-1]) > 1e3 else np.nan,
        "decay_r": float(d["r"]) if d is not None else np.nan,
        "decay_x": float(d["x"]) if d is not None else np.nan,
        "decay_y": float(d["y"]) if d is not None else np.nan,
        "decay_frac": 100.0 * float(d["tau"]) / float(df["tau"].iloc[-1]) if d is not None else np.nan,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csvs", nargs="+", help="two or more CSV files from monograph_plunge")
    p.add_argument("--labels", nargs="*", default=None, help="one label per CSV")
    p.add_argument("--out", default="compare_spins.png", help="output image")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    args = p.parse_args()

    labels = args.labels or [os.path.splitext(os.path.basename(c))[0] for c in args.csvs]
    if len(labels) != len(args.csvs):
        p.error("number of --labels must match number of CSV files")

    apply_style()
    runs = [pd.read_csv(c) for c in args.csvs]
    info = [summarize(df) for df in runs]
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(runs))]

    fig = plt.figure(figsize=(15, 6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.7, 1, 1], wspace=0.28)
    ax = fig.add_subplot(gs[0])

    known = [i["r_plus"] for i in info if not np.isnan(i["r_plus"])]
    r_black = min(known) if known else 1.0
    ax.add_patch(Circle((0, 0), r_black, fc="black", ec=ORANGE, lw=2, zorder=3))
    # Runs with a larger horizon (e.g. a = 0, r+ = 2 M) get a dashed outline in their color.
    for c, inf in zip(colors, info):
        if not np.isnan(inf["r_plus"]) and inf["r_plus"] > r_black + 0.01:
            ax.add_patch(Circle((0, 0), inf["r_plus"], fill=False, ec=c, lw=1.2, ls="--", alpha=0.8, zorder=3))

    for df, lab, c, inf in zip(runs, labels, colors, info):
        ax.plot(df["x"], df["y"], color=c, lw=2.2, label=lab + ("" if np.isnan(inf["r_plus"]) else f"  (horizon r\u208a = {inf['r_plus']:.2f} M)"), zorder=4)
        ax.plot(df["x"].iloc[0], df["y"].iloc[0], "o", color=c, ms=7, zorder=5)
        if not np.isnan(inf["decay_x"]):
            ax.plot(inf["decay_x"], inf["decay_y"], "X", color=c, ms=14, mec="white", mew=1.2, zorder=6)
    lim = 1.1 * max(float(df["r"].max()) for df in runs)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("x (units of M)")
    ax.set_ylabel("y (units of M)")
    ax.set_title("Same launch, different spin  (X = muon decays)")
    ax.legend(loc="lower left", facecolor=BG, edgecolor=GRID, fontsize=9)
    ax.grid(alpha=0.25)

    xs = np.arange(len(runs))
    short = ["\n".join(l.rsplit(" ", 1)) if len(l) > 10 else l for l in labels]

    ax2 = fig.add_subplot(gs[1])
    vals = [i["decay_r"] for i in info]
    bars = ax2.bar(xs, vals, color=colors, width=0.65)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.05, f"{v:.2f}", ha="center", fontsize=10)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(short, fontsize=9)
    ax2.set_ylabel("radius where the muon decays (M)")
    ax2.set_title("Where it dies")
    ax2.set_ylim(0, max(v for v in vals if not np.isnan(v)) * 1.18)

    ax3 = fig.add_subplot(gs[2])
    vals = [i["decay_frac"] for i in info]
    bars = ax3.bar(xs, vals, color=colors, width=0.65)
    for b, v in zip(bars, vals):
        ax3.text(b.get_x() + b.get_width() / 2, v + 1, f"{v:.0f}%", ha="center", fontsize=10)
    ax3.set_xticks(xs)
    ax3.set_xticklabels(short, fontsize=9)
    ax3.set_ylabel("share of the plunge completed at decay (%)")
    ax3.set_title("How far it got")
    ax3.set_ylim(0, 100)

    for a_ in (ax2, ax3):
        for s in ("top", "right"):
            a_.spines[s].set_visible(False)

    fig.savefig(args.out, dpi=160, bbox_inches="tight")

    print(f"{'run':<22}{'tau to horizon':>16}{'t to horizon':>14}{'turns':>8}{'decay r':>10}{'decay %':>9}")
    for lab, inf in zip(labels, info):
        print(f"{lab:<22}{inf['tau_end_us']:>13.3f} us{inf['t_end_us']:>11.3f} us"
              f"{inf['turns']:>8.2f}{inf['decay_r']:>10.3f}{inf['decay_frac']:>8.1f}%")
    print(f"Saved {args.out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
