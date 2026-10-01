"""Reconstruction of fig_network.pdf: network topology at the three
lexicographic anchors (economic / environmental / volume), full instance.

The original script was lost; this rebuilds it from
results/pareto_front/econ_results.json. Fix vs. the lost version: all three
panels share identical xlim/ylim (computed once, from the residential-zone
bounding box) so the equal-aspect scaling is identical across panels -- the
old figure let each panel autoscale independently, which made the left
(economic) panel render a hair larger than the other two.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # raiz del repositorio
import matplotlib.pyplot as plt

from fig_network_common import (
    load_geo, pareto_by_opened, draw_panel, shared_limits, legend_handles,
)

OUT_PATH = ROOT / "figures/fig_network1.pdf"

ANCHORS = [
    ("Economic", [2, 7, 8]),
    ("Environmental", [1, 4, 7, 8, 9, 10, 12, 14]),
    ("Volume", [1, 4, 7, 8, 9, 11, 12, 14]),
]


def main():
    d, zones, hubs, facilities = load_geo()
    pareto = d["pareto"]
    xlim, ylim = shared_limits(zones)

    aspect = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])
    panel_w = 5.2
    fig, axes = plt.subplots(
        1, 3, figsize=(panel_w * 3, panel_w * aspect + 1.3), sharex=True, sharey=True
    )

    for ax, (label, opened) in zip(axes, ANCHORS):
        pt = pareto_by_opened(pareto, opened)
        draw_panel(ax, zones, hubs, facilities, opened, pt["flows_jk"], xlim, ylim)
        ax.set_title(
            f"{label}\n$Z_1$={pt['Z4']:,.0f}€, {pt['n_open']} facilities".replace(",", ","),
            fontsize=11, fontweight="bold",
        )
        ax.set_xlabel("Easting (km, local origin)")

    axes[0].set_ylabel("Northing (km, local origin)")

    fig.suptitle("Network topology at the three mono-objective solutions",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.legend(handles=legend_handles(), loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.06), frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=200, bbox_inches="tight")
    print(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
