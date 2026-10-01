"""Figura 9a: % de recoleccion mensual (colectado / total generado en la REGION)
para las 3 soluciones mono-objetivo (Z4 economica, Z2 ambiental, Z3 cobertura)
mas la solucion de compromiso (punto 21), caso base, un solo panel.

El denominador es la generacion REAL de las 13 458 zonas de la region
(data/instances/zone_full.json -> w_imt), no solo las zonas que quedan conectadas a
la red bajo el D_BAR_IJ del caso base (ver fig9_monthly_collection_sensitivity.py
para la discusion completa de por que esto importa).

Fuentes de datos:
  - data/instances/zone_full.json -> w_imt   (generacion real, TODAS las zonas, fija)
  - results/sensitivity_mono/runs/base/detalle_csv/hubs_temporal.csv  (influjo_kg = colectado
    real hacia hubs; sol 1=Z4, sol 2=Z2, sol 3=Z3, ver .../soluciones.csv)
  - results/sensitivity_compromise/runs/base/detalle_csv/hubs_temporal.csv  (solucion de compromiso, sol=1)
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # raiz del repositorio
import collections
import json

import pandas as pd
import matplotlib.pyplot as plt

ZONE_FULL_JSON = ROOT / "data/instances/zone_full.json"
MONO_DIR = ROOT / "results/sensitivity_mono/runs/base/detalle_csv"
COMP_DIR = ROOT / "results/sensitivity_compromise/runs/base/detalle_csv"
OUT_PATH = ROOT / "figures/fig9_monthly_collection_base.pdf"

MONO_SOL_LABELS = {
    "mono_Z4": ("Mono-Z4 (Economic)", "#c0392b"),
    "mono_Z2": ("Mono-Z2 (Environmental)", "#27ae60"),
    "mono_Z3": ("Mono-Z3 (Coverage)", "#8e44ad"),
}
COMPROMISE_LABEL = ("Compromise (sol. 21)", "#2c5f8a")

plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 300,
})


def true_monthly_generation_kg():
    """Generacion REAL de toda la region (13 458 zonas), por mes."""
    with open(ZONE_FULL_JSON, encoding="utf-8") as fh:
        data = json.load(fh)
    by_month = collections.defaultdict(float)
    for _i, _m, t, w in data["w_imt"]:
        by_month[t] += w
    return pd.Series({t: by_month[t] for t in range(1, 13)}).sort_index()


def monthly_collection_by_sol(detalle_dir):
    df = pd.read_csv(f"{detalle_dir}/hubs_temporal.csv")
    return df.groupby(["sol", "t"])["influjo_kg"].sum()


def mono_sol_map(detalle_dir):
    sol_df = pd.read_csv(f"{detalle_dir}/soluciones.csv")
    return dict(zip(sol_df["estado"], sol_df["sol"]))


def coverage_pct(detalle_dir, sol_number, generation):
    coll = monthly_collection_by_sol(detalle_dir)
    year1 = coll.loc[sol_number].loc[1:12].sort_index()
    return year1 / generation * 100.0, year1


def collect_lines():
    generation = true_monthly_generation_kg()
    lines = {}

    sol_map = mono_sol_map(MONO_DIR)
    for estado, (label, color) in MONO_SOL_LABELS.items():
        pct, abs_kg = coverage_pct(MONO_DIR, sol_map[estado], generation)
        lines[label] = (pct, abs_kg, color)

    pct, abs_kg = coverage_pct(COMP_DIR, 1, generation)
    lines[COMPROMISE_LABEL[0]] = (pct, abs_kg, COMPROMISE_LABEL[1])

    return lines


def main():
    fig, ax = plt.subplots(figsize=(7, 4.8))

    for label, (ys, abs_kg, color) in collect_lines().items():
        annual_t = abs_kg.sum() / 1000.0
        full_label = f"{label}"
        ax.plot(ys.index, ys.values, "-o", color=color, linewidth=1.8,
                markersize=3.5, markeredgecolor="black", label=full_label)

    ax.set_xlabel("Months")
    ax.set_ylabel("% of Ville de Nancy generated waste collected")
    ax.set_xticks(range(1, 13))
    ax.set_title("Monthly collection coverage (base case for a year)")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT_PATH, bbox_inches="tight")
    print(f"Guardado: {OUT_PATH}")


if __name__ == "__main__":
    main()
