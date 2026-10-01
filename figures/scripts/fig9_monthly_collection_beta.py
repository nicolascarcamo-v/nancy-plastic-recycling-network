"""Figura 9-beta: % de recoleccion mensual (colectado / total generado en la
REGION) para las 3 soluciones mono-objetivo (Z4, Z2, Z3) mas la solucion de
compromiso (punto 21), bajo distintos niveles de BETA_MIN (fraccion minima de
la generacion que DEBE colectarse en cada periodo, restriccion beta_min_c).

Un panel por nivel de beta: 0.00 (base, restriccion inactiva), 0.05, 0.10, 0.20.
Misma regla de porcentaje que fig9_monthly_collection_sensitivity.py: el
denominador es la generacion mensual REAL de todas las zonas de
instances/zone_full.json (w_imt), fija e independiente del escenario.

Fuentes de datos (producidas por run_beta_overnight.ps1):
  - results/sensitivity_beta_mono/runs/<run_id>/detalle_csv/{hubs_temporal,soluciones}.csv
  - results/sensitivity_beta_compromise/runs/<run_id>/detalle_csv/hubs_temporal.csv  (sol=1)

Si algun nivel del compromiso resulta infactible (los eps del punto 21 + el piso
beta pueden no ser alcanzables simultaneamente), esa linea se omite y se anota
"infeasible" en el panel — es un resultado, no un error.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # raiz del repositorio
import collections
import json
import os

import pandas as pd
import matplotlib.pyplot as plt

ZONE_FULL_JSON = ROOT / "data/instances/zone_full.json"
MONO_ROOT = ROOT / "results/sensitivity_beta_mono/runs"
COMP_ROOT = ROOT / "results/sensitivity_beta_compromise/runs"
OUT_PATH = ROOT / "figures/fig9_monthly_collection_beta.pdf"

MONO_SOL_LABELS = {
    "mono_Z4": ("Mono-Z4 (Economic)", "#c0392b"),
    "mono_Z2": ("Mono-Z2 (Environmental)", "#27ae60"),
    "mono_Z3": ("Mono-Z3 (Coverage)", "#8e44ad"),
}
COMPROMISE_LABEL = ("Compromise (sol. 21)", "#2c5f8a")

# (run_id, titulo del panel, nivel beta para la linea de referencia)
LEVELS = [
    ("base",       r"$\beta_{\min}=0.00$ (base)", 0.00),
    ("beta_0p05",  r"$\beta_{\min}=0.05$",        0.05),
    ("beta_0p10",  r"$\beta_{\min}=0.10$",        0.10),
    ("beta_0p20",  r"$\beta_{\min}=0.20$",        0.20),
]

plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 300,
})


def true_monthly_generation_kg():
    """Generacion REAL de toda la region, por mes, fija (mismo criterio que
    fig9_monthly_collection_sensitivity.py)."""
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


def coverage_pct_for_run(detalle_dir, sol_number, generation):
    coll = monthly_collection_by_sol(detalle_dir)
    year1 = coll.loc[sol_number].loc[1:12].sort_index()
    return year1 / generation * 100.0


def main():
    generation = true_monthly_generation_kg()

    fig, axes = plt.subplots(1, len(LEVELS), figsize=(15, 4.6), sharey=True)

    for ax, (run_id, title, beta) in zip(axes, LEVELS):
        mono_dir = f"{MONO_ROOT}/{run_id}/detalle_csv"
        comp_dir = f"{COMP_ROOT}/{run_id}/detalle_csv"

        # 3 soluciones mono-objetivo
        if os.path.exists(f"{mono_dir}/soluciones.csv"):
            smap = mono_sol_map(mono_dir)
            for estado, (lbl, color) in MONO_SOL_LABELS.items():
                if estado not in smap:
                    continue
                pct = coverage_pct_for_run(mono_dir, smap[estado], generation)
                ax.plot(pct.index, pct.values, color=color, linewidth=2.0, zorder=3)
        else:
            ax.text(0.5, 0.5, "mono runs\nnot available", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8, color="grey")

        # solucion de compromiso (sol=1); puede ser infactible bajo el piso beta
        try:
            pct = coverage_pct_for_run(comp_dir, 1, generation)
            ax.plot(pct.index, pct.values, color=COMPROMISE_LABEL[1],
                    linewidth=2.0, zorder=3)
        except (FileNotFoundError, KeyError):
            ax.text(0.5, 0.06, "Compromise: infeasible", transform=ax.transAxes,
                    ha="center", va="bottom", fontsize=7,
                    color=COMPROMISE_LABEL[1])

        # piso regulatorio: beta_min * generacion => beta_min en % es una recta
        if beta > 0:
            ax.axhline(beta * 100.0, color="black", linestyle="--",
                       linewidth=1.0, alpha=0.7, zorder=2)
            ax.annotate("floor", (12.0, beta * 100.0), fontsize=7,
                        ha="right", va="bottom", color="black")

        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Month of year")
        ax.set_xticks(range(1, 13))

    axes[0].set_ylabel("% of regionwide\ngenerated waste collected")

    handles = [plt.Line2D([0], [0], color=c, lw=2, label=lbl)
               for lbl, c in list(MONO_SOL_LABELS.values()) + [COMPROMISE_LABEL]]
    handles.append(plt.Line2D([0], [0], color="black", lw=1, linestyle="--",
                              label=r"Regulatory floor $\beta_{\min}$"))
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, -0.06), fontsize=8.5)

    fig.suptitle("Monthly collection coverage vs. true regionwide waste generation\n"
                 r"under increasing minimum-collection floors ($\beta_{\min}$)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PATH, bbox_inches="tight")
    print(f"Guardado: {OUT_PATH}")


if __name__ == "__main__":
    main()
