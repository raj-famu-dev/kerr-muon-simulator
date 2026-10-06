#!/usr/bin/env python3
"""Coordinate time t versus proper time tau: the divergence at the horizon.

Panels:
  1. t against tau (the muon's clock against the distant observer's clock)
  2. dt/dtau against tau, on a log axis: the rate blows up at the horizon
  3. (only with --spin) t against -ln(r - r+): a straight line means t diverges
     logarithmically, with the slope (r+^2 + a^2)/(r+ - r-) predicted by theory

Usage:
    python plot_time_dilation.py plunge.csv --spin 0.9 --out time_dilation.png
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
CYAN, ORANGE, RED, YELLOW = "#6ee7ff", "#ff9d2e", "#ff4d6d", "#ffd27a"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
        "axes.edgecolor": GRID, "axes.labelcolor": FG, "xtick.color": FG,
        "ytick.color": FG, "text.color": FG, "grid.color": GRID, "font.size": 11,
    })


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="CSV written by monograph_plunge")
    p.add_argument("--spin", type=float, default=None,
                   help="black hole spin a; adds the logarithmic-divergence panel")
    p.add_argument("--out", default=None, help="output image (default: next to the CSV)")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    args = p.parse_args()

    apply_style()
    df = pd.read_csv(args.csv)
    tau, t = df["tau_us"].to_numpy(), df["t_us"].to_numpy()
    dead = df[df["muon_alive"] == 0]
    tau_decay = float(dead.iloc[0]["tau_us"]) if len(dead) else None

    npan = 3 if args.spin is not None else 2
    fig, axes = plt.subplots(1, npan, figsize=(6.2 * npan, 5.2))

    # ---- panel 1: t against tau ------------------------------------------
    ax = axes[0]
    ax.plot(tau, tau, "--", color="#aab", lw=1.2, label="t = \u03c4  (no gravity)")
    ax.plot(tau, t, color=ORANGE, lw=3, label="distant observer's clock t")
    if tau_decay is not None:
        ax.axvline(tau_decay, color=YELLOW, ls=":", lw=1.5)
        ax.text(tau_decay, 0.62 * t.max(), " muon decays\n (2.2 \u03bcs)", color=YELLOW,
                va="top", ha="left", fontsize=10)
        ax.axvspan(tau_decay, tau.max(), color="white", alpha=0.04)
    ax.set_xlabel("muon proper time \u03c4 (\u03bcs)")
    ax.set_ylabel("coordinate time t (\u03bcs)")
    ax.set_title("t races ahead of \u03c4")
    ax.legend(loc="upper left", facecolor=BG, edgecolor=GRID, fontsize=9)
    ax.grid(alpha=0.25)
    ax.text(0.04, 0.58, f"at the horizon:\n\u03c4 = {tau[-1]:.2f} \u03bcs\nt = {t[-1]:.2f} \u03bcs\n"
                        f"({t[-1] / tau[-1]:.1f}\u00d7 more,\nstill growing)",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color=FG)

    # ---- panel 2: dt/dtau ---------------------------------------------------
    ax = axes[1]
    ax.semilogy(tau, df["dt_dtau"].to_numpy(), color=CYAN, lw=2.5)
    if tau_decay is not None:
        ax.axvline(tau_decay, color=YELLOW, ls=":", lw=1.5)
    ax.set_xlabel("muon proper time \u03c4 (\u03bcs)")
    ax.set_ylabel("dt / d\u03c4  (distant seconds per muon second)")
    ax.set_title("the clock-rate ratio blows up")
    ax.grid(alpha=0.25, which="both")

    # ---- panel 3: logarithmic divergence -------------------------------------
    if args.spin is not None:
        a = args.spin
        r_plus = 1 + np.sqrt(1 - a * a)
        r_minus = 1 - np.sqrt(1 - a * a)
        gap = df["r"].to_numpy() - r_plus
        keep = gap > 0
        xlog = -np.log(gap[keep])
        tm = df["t"].to_numpy()[keep]
        ax = axes[2]
        ax.plot(xlog, tm, color=ORANGE, lw=3, label="engine")
        slope = (r_plus ** 2 + a * a) / (r_plus - r_minus)
        i = int(0.9 * len(xlog))
        xs = np.linspace(xlog[i], xlog[-1], 2)
        ax.plot(xs, tm[i] + slope * (xs - xlog[i]), "--", color="white", lw=1.5,
                label=f"theory: slope = (r\u208a\u00b2+a\u00b2)/(r\u208a\u2212r\u208b) = {slope:.2f}")
        ax.set_xlabel("\u2212ln(r \u2212 r\u208a)   (grows as the muon nears the horizon)")
        ax.set_ylabel("coordinate time t (units of M)")
        ax.set_title("t diverges logarithmically")
        ax.legend(loc="upper left", facecolor=BG, edgecolor=GRID, fontsize=9)
        ax.grid(alpha=0.25)

    out = args.out or str(args.csv).rsplit(".", 1)[0] + "_time_dilation.png"
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print(f"Saved {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
