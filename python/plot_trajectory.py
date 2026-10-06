#!/usr/bin/env python3
"""Orbital path of the muon, colored by its own proper time.

The small inset magnifies the last turns around the horizon, where frame
dragging winds the path around the black hole.

Usage:
    python plot_trajectory.py plunge.csv --out trajectory.png
"""
import argparse
import sys

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Circle

BG, FG, GRID = "#0b0720", "#e8e0ff", "#2a2250"
CYAN, ORANGE, RED, PURPLE = "#6ee7ff", "#ff9d2e", "#ff4d6d", "#c04bff"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
        "axes.edgecolor": GRID, "axes.labelcolor": FG, "xtick.color": FG,
        "ytick.color": FG, "text.color": FG, "grid.color": GRID, "font.size": 11,
    })


def find_horizon(df, spin=None):
    """Horizon radius: exact if the spin is given; otherwise read from the data when
    the run ended at the horizon (dt/dtau has blown up). Returns None if unknown."""
    if spin is not None:
        return 1.0 + float(np.sqrt(1.0 - spin * spin))
    if float(df["dt_dtau"].iloc[-1]) > 1e3:
        return float(df["r"].iloc[-1])
    return None


def bright_plasma():
    """Plasma without its near-black start, so the early path stays visible."""
    return ListedColormap(plt.cm.plasma(np.linspace(0.2, 1.0, 256)))


def draw_hole(ax, r_plus, label=True):
    if r_plus is None:
        return
    if r_plus < 1.99:  # the ergosphere (r = 2 M on the equator) exists only if the hole spins
        ax.add_patch(Circle((0, 0), 2.0, fill=False, ec=PURPLE, lw=1.2, ls="--", alpha=0.9, zorder=2))
        if label:
            ax.text(0, 2.0 + 0.12, "ergosphere", color=PURPLE, ha="center", fontsize=9)
    ax.add_patch(Circle((0, 0), r_plus, fc="black", ec=ORANGE, lw=2.2, zorder=3))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="CSV written by monograph_plunge")
    p.add_argument("--out", default=None, help="output image (default: next to the CSV)")
    p.add_argument("--title", default=None)
    p.add_argument("--spin", type=float, default=None,
                   help="black hole spin a (needed only if the run did not reach the horizon)")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    args = p.parse_args()

    apply_style()
    df = pd.read_csv(args.csv)
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    tau_us = df["tau_us"].to_numpy()
    r_plus = find_horizon(df, args.spin)
    if r_plus is None:
        print("Note: the run did not end at the horizon; pass --spin to draw it.")

    pts = np.column_stack([x, y]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    cmap, norm = bright_plasma(), Normalize(tau_us.min(), tau_us.max())

    def make_lc(width):
        lc = LineCollection(segs, cmap=cmap, norm=norm, linewidth=width, zorder=4)
        lc.set_array(tau_us[:-1])
        return lc

    fig, ax = plt.subplots(figsize=(8.5, 8))
    draw_hole(ax, r_plus)
    lc = make_lc(2.6)
    ax.add_collection(lc)

    ax.plot(x[0], y[0], "o", color=CYAN, ms=9, zorder=6)
    ax.annotate("launch", (x[0], y[0]), xytext=(0, -22), textcoords="offset points",
                color=CYAN, ha="center")

    dead = df[df["muon_alive"] == 0]
    if len(dead):
        d = dead.iloc[0]
        ax.plot(d["x"], d["y"], "X", color=RED, ms=15, mec="white", mew=1.2, zorder=7)
        ax.annotate(f"muon decays\n\u03c4 = {d['tau_us']:.2f} \u03bcs, r = {d['r']:.2f} M",
                    (d["x"], d["y"]), xytext=(0, 66), textcoords="offset points", color=RED,
                    ha="center", arrowprops=dict(arrowstyle="->", color=RED), fontsize=10)

    lim = 1.12 * float(df["r"].max())
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("x (units of M)")
    ax.set_ylabel("y (units of M)")
    ax.grid(alpha=0.25)

    # ---- inset: how the orbital angle accumulates with proper time ----------
    axin = ax.inset_axes([0.64, 0.08, 0.32, 0.24])
    axin.set_facecolor(BG)
    axin.plot(tau_us, df["phi"].to_numpy() / (2 * np.pi), color=ORANGE, lw=2)
    axin.set_xlabel("muon proper time \u03c4 (\u03bcs)", fontsize=8, labelpad=2)
    axin.set_ylabel("angle swept (turns)", fontsize=8, labelpad=2)
    axin.tick_params(labelsize=8, length=2)
    axin.grid(alpha=0.25)
    for sp in axin.spines.values():
        sp.set_edgecolor(PURPLE)
    axin.set_title("frame dragging winds it up", fontsize=9, color=PURPLE, pad=3)

    cb = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("muon proper time \u03c4 (\u03bcs)")
    cb.outline.set_edgecolor(GRID)

    title = args.title or "Muon plunge into a rotating black hole"
    ax.set_title(title, pad=24, fontsize=14)
    ax.text(0.5, 1.02, f"muon clock: {df['tau_us'].iloc[-1]:.2f} \u03bcs   |   distant clock: "
                       f"{df['t_us'].iloc[-1]:.2f} \u03bcs" + ("   (until it reaches the horizon)" if float(df["dt_dtau"].iloc[-1]) > 1e3 else ""),
            transform=ax.transAxes, ha="center", fontsize=10, color="#b9b3d9")

    out = args.out or str(args.csv).rsplit(".", 1)[0] + "_trajectory.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print(f"Saved {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
