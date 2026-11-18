"""Plot Fig. 8 (mutual information) and Figs. 13-14 (BLER).

Run:
    python plots.py mi --mode smoke     # tiny check, writes results/fig8_mutual_information.png
    python plots.py mi --mode quick     # rough curves (tens of minutes to hours)
    python plots.py bler --mode smoke   # one BLER point, writes results/fig13_bler.png and fig14
    python plots.py bler --mode quick
    python plots.py all --mode quick
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from experiments import get_cfg, mi_curve, bler_curve


COLORS = [
    (0.58, 0.00, 0.83),
    (0.00, 0.62, 0.45),
    (0.34, 0.71, 0.91),
    (0.90, 0.62, 0.00),
    (0.94, 0.89, 0.26),
    (0.00, 0.45, 0.70),
    (0.84, 0.37, 0.00),
    (0.00, 0.00, 0.00),
]
MARKERS = ["+", "x", "*", "s", "o", "d", "^", "v"]
TITLES = [
    "CH-1: DNA channel model",
    "CH-2: Deletion channel",
    "CH-3: Insertion channel",
    "CH-4: Substitution channel",
    "CH-5: IDS channel",
]


def _style(k):
    return COLORS[k % len(COLORS)], MARKERS[k % len(MARKERS)]


def _save(fig, cfg, name):
    path = Path(cfg.cache_dir) / f"{name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {path}")


def fig8_mutual_information(cfg):
    """Paper Fig. 8: I(C;U) vs alpha for CH-1..5, alphabet w41."""
    fig, axes = plt.subplots(3, 2, figsize=(8.5, 10))
    axes = axes.ravel()
    for ch in range(1, 6):
        ax = axes[ch - 1]
        for k, r in enumerate(cfg.mi.r_list):
            color, marker = _style(k)
            I = mi_curve(cfg, ch, "w41", r, cfg.mi.alphas)
            ax.plot(cfg.mi.alphas, I, "-" + marker, color=color, markersize=5, label=f"r={r}")
        ax.set_xlim(0, 1.4)
        ax.set_ylim(0, 7)
        ax.grid(True)
        ax.set_title(TITLES[ch - 1])
        ax.set_xlabel("Error scale alpha")
        ax.set_ylabel("Mutual information (bit)")
    axes[5].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[5].legend(handles, labels, loc="center left", title="sequencing depth")
    fig.tight_layout()
    _save(fig, cfg, "fig8_mutual_information")


def _plot_bler(ax, S, color, linestyle):
    ok = np.isfinite(S["nerr"])
    if not np.any(ok):
        return
    b = S["nerr"][ok] / S["nblk"][ok]
    a = S["alphas"][ok]
    pos = b > 0
    if np.any(pos):
        ax.plot(a[pos], b[pos], linestyle, color=color, marker="o", markersize=4)
    if np.any(~pos):
        ax.plot(a[~pos], 1.0 / S["nblk"][ok][~pos], linestyle="none",
                color=color, marker="v", markersize=5)


def fig13_14_bler(cfg):
    """Paper Figs. 13-14: BLER of the concatenated scheme (CH-1, w41)."""
    fig, ax = plt.subplots(figsize=(7, 5.2))
    for i, M in enumerate(cfg.bler.Ms):
        color, _ = _style(i)
        for r in cfg.bler.r_list:
            S = bler_curve(cfg, M, r, cfg.bler.alphas)
            ls = "-" if r == max(cfg.bler.r_list) else "--"
            _plot_bler(ax, S, color, ls)
            ax.plot([], [], ls, color=color, label=f"M={M}, r={r}")
    ax.set_yscale("log")
    ax.set_xlim(0.2, 1.4)
    ax.set_ylim(1e-5, 1)
    ax.grid(True, which="both")
    ax.set_xlabel("Error scale alpha")
    ax.set_ylabel("Block error rate (BLER)")
    ax.set_title("BLER vs error scale (solid: largest r, dashed: smaller r)")
    ax.legend(loc="lower right")
    fig.tight_layout()
    _save(fig, cfg, "fig13_bler")

    fig, ax = plt.subplots(figsize=(7, 5.2))
    for j, r in enumerate(cfg.bler.r14_list):
        color, _ = _style(j)
        S = bler_curve(cfg, 768, r, cfg.bler.alphas)
        _plot_bler(ax, S, color, "-")
        ax.plot([], [], "-", color=color, label=f"r={r}")
    ax.set_yscale("log")
    ax.set_xlim(0.2, 1.4)
    ax.set_ylim(1e-5, 1)
    ax.grid(True, which="both")
    ax.set_xlabel("Error scale alpha")
    ax.set_ylabel("Block error rate (BLER)")
    ax.set_title("BLER vs error scale, M = 768")
    ax.legend(loc="lower right")
    fig.tight_layout()
    _save(fig, cfg, "fig14_bler_vs_r")


def main():
    parser = argparse.ArgumentParser(description="Plot mutual information and BLER curves")
    parser.add_argument("which", choices=["mi", "bler", "all"], help="which figures to run")
    parser.add_argument("--mode", default="smoke", choices=["smoke", "quick", "full"])
    args = parser.parse_args()
    cfg = get_cfg(args.mode)
    print(f"mode={args.mode}  n={cfg.n}  cache={cfg.cache_dir}")
    if args.which in ("mi", "all"):
        print("== Fig. 8 mutual information (w41) ==")
        fig8_mutual_information(cfg)
    if args.which in ("bler", "all"):
        print("== Figs. 13-14 BLER (w41, CH-1) ==")
        fig13_14_bler(cfg)


if __name__ == "__main__":
    main()
