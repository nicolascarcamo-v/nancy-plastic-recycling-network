"""Reconstruction of fig_network_comparison.pdf: network growth from the
strongest reliable mono-objective economic configuration (sol 4) to the
compromise solution (sol 21), full instance.

Original script lost; rebuilt from results/pareto_front/econ_results.json
(facility/hub geometry + flows_jk) and
results/pareto_front/detalle_csv/soluciones.csv (sol metadata: Z4, n_plantas,
gap_pct). Both panels share identical xlim/ylim so panel (a) does not render
larger than panel (b), which was the scaling bug in the lost version.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # raiz del repositorio
import matplotlib.pyplot as plt

from fig_network_common import (
    load_geo, pareto_by_opened, draw_panel, shared_limits, legend_handles,
    read_soluciones_csv, sol_by_id,
)

OUT_PATH = ROOT / "figures/fig_network_comparison1.pdf"

PANELS = [
    ("Best Reliable (gap-wise) Economic Solution", 4, [1, 4, 7, 8]),
    ("Compromise Solution", 21, [1, 4, 7, 8, 9, 10, 14]),
]


def main():
    d, zones, hubs, facilities = load_geo()
    pareto = d["pareto"]
    rows = read_soluciones_csv()
    xlim, ylim = shared_limits(zones)

    aspect = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])
    panel_w = 6.4
    fig, axes = plt.subplots(
        1, 2, figsize=(panel_w * 2, panel_w * aspect + 1.3), sharex=True, sharey=True
    )

    for ax, (label, sol_id, opened) in zip(axes, PANELS):
        row = sol_by_id(rows, sol_id)
        pt = pareto_by_opened(pareto, opened)
        draw_panel(ax, zones, hubs, facilities, opened, pt["flows_jk"], xlim, ylim)
        z1 = float(row["Z4_NPV_EUR"])
        gap = float(row["gap_pct"])
        ax.set_title(
            f"{label}\n$Z_1$={z1:,.0f}€, {row['n_plantas']} facilities, gap={gap:.1f}%",
            fontsize=11, fontweight="bold",
        )
        ax.set_xlabel("Easting (km)")

    axes[0].set_ylabel("Northing (km)")

    fig.suptitle(
        "Strongest (and reliable) Economic Solution vs Compromise Solution",
        fontsize=13, fontweight="bold", y=1.03,
    )
    fig.legend(handles=legend_handles(), loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.06), frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=200, bbox_inches="tight")
    print(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
