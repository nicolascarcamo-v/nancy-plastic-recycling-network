"""Figura 9b: replica exacta del estilo de fig9_monthly_collection_base.py (lineas
con marcadores, leyenda limpia) para cada uno de los 4 escenarios de
sensibilidad OFAT, mostrando un panel POR NIVEL (bajo / base / alto) para no
mezclar niveles dentro de una misma linea.

Genera un PDF por escenario, cada uno con 3 paneles lado a lado (mismo eje Y
para poder comparar magnitudes entre niveles a simple vista):
  fig9_monthly_collection_price.pdf
  fig9_monthly_collection_catchment.pdf
  fig9_monthly_collection_muncredit.pdf
  fig9_monthly_collection_facility.pdf

catchment_0p10 no tiene solucion de compromiso factible (ver
results/sensitivity_compromise/sensibilidad_compromiso.csv, factible=0); ese panel simplemente
muestra las 3 soluciones mono-obj y una nota indicando la inviabilidad.

Fuentes de datos: identicas a fig9_monthly_collection_base.py, una carpeta de
corrida (run_id) por nivel en vez de siempre "base".
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # raiz del repositorio
import collections
import json

import pandas as pd
import matplotlib.pyplot as plt

ZONE_FULL_JSON = ROOT / "data/instances/zone_full.json"
MONO_ROOT = ROOT / "results/sensitivity_mono/runs"
COMP_ROOT = ROOT / "results/sensitivity_compromise/runs"

MONO_SOL_LABELS = {
    "mono_Z4": ("Mono-Z1 (Economic)", "#c0392b"),
    "mono_Z2": ("Mono-Z2 (Environmental)", "#27ae60"),
    "mono_Z3": ("Mono-Z3 (Coverage)", "#8e44ad"),
}
COMPROMISE_LABEL = ("Compromise (sol. 21)", "#2c5f8a")

SCENARIOS = [
    ("price", ROOT / "figures/fig9_monthly_collection_price.pdf", "Price volatility",
     [("Low (0.7x)", "price_0p70"), ("Base (1.0x)", "base"), ("High (1.3x)", "price_1p30")]),
    ("catchment", ROOT / "figures/fig9_monthly_collection_catchment.pdf", "Catchment sensitivity",
     [("Low (0.10 km)", "catchment_0p10"), ("Base (0.25 km)", "base"), ("High (0.50 km)", "catchment_0p50")]),
    ("muncredit", ROOT / "figures/fig9_monthly_collection_muncredit.pdf", "Municipal credit intensity",
     [("Low (0.00)", "muncredit_0p00"), ("Base (0.25)", "base"), ("High (0.50)", "muncredit_0p50")]),
    ("facility", ROOT / "figures/fig9_monthly_collection_facility.pdf", "Facility cost/capacity",
     [("Low (0.7x)", "facility_0p70"), ("Base (1.0x)", "base"), ("High (1.3x)", "facility_1p30")]),
]

plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 300,
})


def true_monthly_generation_kg():
    """Generacion REAL de toda la region (13 458 zonas), por mes; fija e
    independiente del escenario/nivel."""
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
    return year1 / generation * 100.0


def collect_lines_for_run(run_id, generation):
    """Return dict {label: (pct_series, color)} for one run_id, in the same
    style as fig9_monthly_collection_base.collect_lines(). Skips the
    compromise line if that run has no feasible detalle_csv."""
    lines = {}
    mono_dir = f"{MONO_ROOT}/{run_id}/detalle_csv"
    sol_map = mono_sol_map(mono_dir)
    for estado, (label, color) in MONO_SOL_LABELS.items():
        pct = coverage_pct(mono_dir, sol_map[estado], generation)
        lines[label] = (pct, color)

    comp_dir = f"{COMP_ROOT}/{run_id}/detalle_csv"
    try:
        pct = coverage_pct(comp_dir, 1, generation)
        lines[COMPROMISE_LABEL[0]] = (pct, COMPROMISE_LABEL[1])
    except FileNotFoundError:
        lines["_compromise_infeasible"] = True

    return lines


def plot_panel(ax, lines, panel_title):
    infeasible = lines.pop("_compromise_infeasible", False)
    for label, (ys, color) in lines.items():
        ax.plot(ys.index, ys.values, "-o", color=color, linewidth=1.8,
                markersize=3.5, markeredgecolor="black", label=label)

    if infeasible:
        ax.text(0.5, 0.06, "Compromise: infeasible", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=7.5, color=COMPROMISE_LABEL[1])

    ax.set_xlabel("Months")
    ax.set_xticks(range(1, 13))
    ax.set_title(panel_title, fontsize=9.5)


def make_combined_figure(generation, out_path=ROOT / "figures/fig9_monthly_collection_all_scenarios.pdf"):
    """Una sola figura: 4 filas (escenarios) x 3 columnas (bajo/base/alto),
    mismo eje Y en toda la figura para comparar cualquier panel entre si."""
    n_rows = len(SCENARIOS)
    fig, axes = plt.subplots(n_rows, 3, figsize=(13, 3.6 * n_rows), sharey=True, sharex=True)

    legend_handles, legend_labels = None, None
    for row, (sid, _out_path, scenario_title, level_runs) in enumerate(SCENARIOS):
        for col, (level_label, run_id) in enumerate(level_runs):
            ax = axes[row, col]
            lines = collect_lines_for_run(run_id, generation)
            plot_panel(ax, lines, level_label)
            ax.set_xlabel("Months" if row == n_rows - 1 else "")
            if row < n_rows - 1:
                ax.tick_params(labelbottom=False)
            h, l = ax.get_legend_handles_labels()
            if legend_handles is None or len(h) > len(legend_handles):
                legend_handles, legend_labels = h, l

        axes[row, 0].set_ylabel(f"{scenario_title}\n% collected", fontsize=8.5)

    fig.legend(legend_handles, legend_labels, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.01), fontsize=9)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.35)
    fig.savefig(out_path, bbox_inches="tight")
    print(f"Guardado: {out_path}")


def main():
    generation = true_monthly_generation_kg()

    for sid, out_path, scenario_title, level_runs in SCENARIOS:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)

        for ax, (level_label, run_id) in zip(axes, level_runs):
            lines = collect_lines_for_run(run_id, generation)
            plot_panel(ax, lines, level_label)

        axes[0].set_ylabel("% of Ville de Nancy generated waste collected")

        handles, labels = axes[1].get_legend_handles_labels()
        if len(handles) < 4:
            for ax in axes:
                h, l = ax.get_legend_handles_labels()
                if len(h) > len(handles):
                    handles, labels = h, l
        fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
                   bbox_to_anchor=(0.5, -0.05), fontsize=8.5)

        fig.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        print(f"Guardado: {out_path}")

    make_combined_figure(generation)


if __name__ == "__main__":
    main()
