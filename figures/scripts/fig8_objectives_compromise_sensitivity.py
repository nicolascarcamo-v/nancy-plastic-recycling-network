"""Figura 8: cómo se mueven los 4 objetivos (Z1..Z4) de la solución de
compromiso (punto 21) bajo el análisis de sensibilidad OFAT.

Fuente de datos: results/sensitivity_compromise/sensibilidad_compromiso.csv, generado por
sensibilidad_compromiso.py. Cada fila es una re-optimización del subproblema
AUGMECON del punto de compromiso con un parámetro perturbado, ceteris paribus.

Para poder superponer los 4 objetivos (EUR, kgCO2e, % cobertura, EUR) en un
mismo eje, cada uno se expresa como variación porcentual respecto a su propio
valor base (fila con es_base=1 dentro de cada escenario).
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # raiz del repositorio
import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = ROOT / "results/sensitivity_compromise/sensibilidad_compromiso.csv"
OUT_PATH = ROOT / "figures/fig8_objectives_compromise_sensitivity.pdf"

OBJ_COLS = {
    "Z1_VAN_peor_EUR": ("Z1 (VAN)", "#c0392b"),
    "Z2_GHG_kgCO2e":   ("Z2 (Environmenta)",            "#27ae60"),
    "Z3_cobertura":    ("Z3 (Coberture)",      "#8e44ad"),
    "Z4_NPV_EUR":      ("Compromise solution", "#2c5f8a"),
}

SCENARIOS = [
    ("price",     "Price volatility",           "PRICE_FACTOR"),
    ("catchment", "Catchment sensitivity",       "D_BAR_IJ"),
    ("muncredit", "Municipal credit intensity",  "MUNICIPAL_CREDIT_SHARE"),
    ("facility",  "Facility cost/capacity",      "FACILITY_FACTOR"),
]

plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 300,
})


def pct_change(series, base_value):
    return (series - base_value) / abs(base_value) * 100.0


def main():
    df = pd.read_csv(CSV_PATH)

    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), sharey=True)

    for ax, (sid, label, param) in zip(axes, SCENARIOS):
        sub = df[df["escenario"] == sid].sort_values("nivel", na_position="first")
        base_row = sub[sub["es_base"] == 1].iloc[0]

        feasible = sub[sub["factible"] == 1]
        infeasible_levels = sub.loc[sub["factible"] == 0, "nivel"].tolist()

        for col, (obj_label, color) in OBJ_COLS.items():
            base_val = base_row[col]
            xs = feasible["nivel"].tolist()
            ys = pct_change(feasible[col], base_val).tolist()
            ax.plot(xs, ys, "o-", color=color, markeredgecolor="black",
                    markersize=4, linewidth=1.3, label=obj_label)

        ax.axhline(0, ls=":", color="gray", lw=1)

        for lv in infeasible_levels:
            ax.axvline(lv, color="firebrick", ls="--", lw=1)
            ax.text(lv, ax.get_ylim()[0], "infeasible", rotation=90,
                     color="firebrick", fontsize=7, va="bottom", ha="right")

        ax.set_title(label, fontsize=8.5)
        ax.set_xlabel(param, fontsize=8)

    axes[0].set_ylabel(r"$\Delta$ objective vs. base (%)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.06), fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PATH, bbox_inches="tight")
    print(f"Guardado: {OUT_PATH}")


if __name__ == "__main__":
    main()
