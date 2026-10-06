#!/usr/bin/env python3
"""Animated GIF of the plunge, with the muon's clock and the distant clock.

Frames are evenly spaced in the muon's proper time. The distant observer's
clock t is read from the same data, so you can watch it pull ahead. The trail
is drawn from the full-resolution data, so the spiral near the horizon stays smooth.

Usage:
    python animate_plunge.py plunge.csv --out plunge.gif
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Circle

BG, FG, GRID = "#0b0720", "#e8e0ff", "#2a2250"
CYAN, ORANGE, RED, PURPLE = "#6ee7ff", "#ff9d2e", "#ff4d6d", "#c04bff"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
        "axes.edgecolor": GRID, "axes.labelcolor": FG, "xtick.color": FG,
        "ytick.color": FG, "text.color": FG, "font.size": 11,
    })


def find_horizon(df, spin=None):
    """Horizon radius: exact if the spin is given; otherwise read from the data when
    the run ended at the horizon (dt/dtau has blown up). Returns None if unknown."""
    if spin is not None:
        return 1.0 + float(np.sqrt(1.0 - spin * spin))
    if float(df["dt_dtau"].iloc[-1]) > 1e3:
        return float(df["r"].iloc[-1])
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="CSV written by monograph_plunge")
    p.add_argument("--out", default="plunge.gif", help="output GIF")
    p.add_argument("--frames", type=int, default=150, help="number of animation frames")
    p.add_argument("--fps", type=int, default=25)
    p.add_argument("--hold", type=int, default=25, help="extra frames held at the end")
    p.add_argument("--title", default="Muon falling into a Kerr black hole")
    p.add_argument("--spin", type=float, default=None,
                   help="black hole spin a (needed only if the run did not reach the horizon)")
    args = p.parse_args()

    apply_style()
    df = pd.read_csv(args.csv)
    tau = df["tau"].to_numpy()
    x_f, y_f = df["x"].to_numpy(), df["y"].to_numpy()
    tau_us_f = df["tau_us"].to_numpy()

    # Frames evenly spaced in proper time.
    tau_u = np.linspace(0.0, tau[-1], args.frames)
    idx = np.clip(np.searchsorted(tau, tau_u), 1, len(tau) - 1)
    x, y = np.interp(tau_u, tau, x_f), np.interp(tau_u, tau, y_f)
    tau_us = np.interp(tau_u, tau, tau_us_f)
    t_us = np.interp(tau_u, tau, df["t_us"].to_numpy())

    r_plus = find_horizon(df, args.spin)
    dead = df[df["muon_alive"] == 0]
    tau_decay_us = float(dead.iloc[0]["tau_us"]) if len(dead) else np.inf

    pts = np.column_stack([x_f, y_f]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)

    # View fitted to the path, so the animation is not mostly empty space.
    pad = 0.7
    x0, x1 = min(x_f.min(), -2.2) - pad, max(x_f.max(), 2.2) + pad
    y0, y1 = min(y_f.min(), -2.2) - pad, max(y_f.max(), 2.2) + pad
    fig_w, ax_w = 7.0, 6.7
    ax_h = ax_w * (y1 - y0) / (x1 - x0)
    head = 1.35
    fig_h = ax_h + head
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0.02, 0.02 / fig_h * 1.0, 0.96, ax_h / fig_h])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.axis("off")

    if r_plus is not None:
        if r_plus < 1.99:
            ax.add_patch(Circle((0, 0), 2.0, fill=False, ec=PURPLE, lw=1, ls="--", alpha=0.8))
            ax.text(0, 2.12, "ergosphere", color=PURPLE, ha="center", fontsize=8)
        ax.add_patch(Circle((0, 0), r_plus, fc="black", ec=ORANGE, lw=2, zorder=3))

    cmap = ListedColormap(plt.cm.plasma(np.linspace(0.2, 1.0, 256)))
    trail = LineCollection([], cmap=cmap, norm=Normalize(0, tau_us_f[-1]), linewidth=2.4, zorder=4)
    ax.add_collection(trail)
    dot, = ax.plot([], [], "o", color=CYAN, ms=9, zorder=6)
    ring, = ax.plot([], [], "o", mfc="none", mec=CYAN, ms=17, mew=1.2, alpha=0.6, zorder=6)

    def fy(inches_from_top):
        return 1.0 - inches_from_top / fig_h

    fig.text(0.5, fy(0.35), args.title, ha="center", fontsize=14, color=FG, weight="bold")
    txt_tau = fig.text(0.05, fy(0.80), "", fontsize=12, color=CYAN, family="monospace")
    txt_t = fig.text(0.05, fy(1.12), "", fontsize=12, color=ORANGE, family="monospace")
    txt_state = fig.text(0.95, fy(0.80), "", fontsize=12, color=CYAN, ha="right", weight="bold")

    total = args.frames + args.hold

    def update(i):
        k = min(i, args.frames - 1)
        n = int(idx[k])
        trail.set_segments(segs[:n])
        trail.set_array(tau_us_f[:n])
        alive = tau_us[k] < tau_decay_us - 1e-9
        col = CYAN if alive else RED
        dot.set_data([x[k]], [y[k]])
        dot.set_color(col)
        ring.set_data([x[k]], [y[k]])
        ring.set_markeredgecolor(col)
        txt_tau.set_text(f"muon clock    \u03c4 = {tau_us[k]:6.2f} \u03bcs")
        txt_t.set_text(f"distant clock t = {t_us[k]:6.2f} \u03bcs")
        txt_state.set_text("ALIVE" if alive else "DECAYED")
        txt_state.set_color(col)
        return trail, dot, ring, txt_tau, txt_t, txt_state

    ani = FuncAnimation(fig, update, frames=total, blit=False)
    ani.save(args.out, writer=PillowWriter(fps=args.fps), dpi=80)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
