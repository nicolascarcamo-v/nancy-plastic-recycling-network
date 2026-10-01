"""Figura 9: % de recoleccion mensual (colectado / total generado en la REGION)
para las 3 soluciones mono-objetivo (Z4 economica, Z2 ambiental, Z3 cobertura)
mas la solucion de compromiso (punto 21), un panel por escenario de sensibilidad
OFAT (price, catchment, muncredit, facility), mostrando los 3 niveles de cada
escenario (bajo / base / alto).

*** CORRECCION IMPORTANTE DE LA REGLA DE PORCENTAJE ***
La columna oficial `cobertura_pct` de sensibilidad_monoobj.csv /
sensibilidad_compromiso.csv se calcula como colectado / generado, pero
"generado" ahi es la suma de `generado_kg` en zonas_cobertura.csv, que solo
incluye las zonas que quedaron CONECTADAS a la red bajo el D_BAR_IJ (radio de
captacion) de ESE escenario/nivel (ver sensibilidad.py:_filter_connected_graph,
que descarta del grafo cualquier zona i sin arco a <= D_BAR_IJ). Al achicar el
radio, miles de zonas quedan fuera del modelo por completo y su generacion
DESAPARECE del denominador, inflando artificialmente el % de cobertura.
Ejemplo real: catchment_0p10 (D_BAR_IJ=0.10 km) deja solo 386 de las 13 458
zonas de la region conectadas; el modelo "cubre" el 98% de esas 386 zonas,
pero eso es apenas el 3.2% de la generacion REAL de toda la region.

Por eso aqui el denominador es fijo e independiente del escenario/nivel: la
generacion mensual de TODAS las 13 458 zonas de data/instances/zone_full.json
(w_imt), tal cual entra al modelo antes de cualquier filtro de conectividad.
Esto es lo correcto para dar recomendaciones regulatorias: muestra cuanto de
la basura que REALMENTE se genera en la region queda fuera de la cadena
formal (y por lo tanto es candidata a reciclaje paralelo/informal).

Fuentes de datos:
  - data/instances/zone_full.json -> w_imt   (generacion real, TODAS las zonas, fija)
  - results/sensitivity_mono/runs/<run_id>/detalle_csv/hubs_temporal.csv  (influjo_kg = colectado
    real hacia hubs; sol 1=Z4, sol 2=Z2, sol 3=Z3, ver .../soluciones.csv)
  - results/sensitivity_compromise/runs/<run_id>/detalle_csv/hubs_temporal.csv  (solucion de compromiso, sol=1)

catchment_0p10 no tiene solucion de compromiso factible (ver
results/sensitivity_compromise/sensibilidad_compromiso.csv, factible=0) -> esa linea se omite y se
anota "infeasible" en el panel.
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
OUT_PATH = ROOT / "figures/fig9_monthly_collection_sensitivity.pdf"

MONO_SOL_LABELS = {
    "mono_Z4": ("Mono-Z4 (Economic)", "#c0392b"),
    "mono_Z2": ("Mono-Z2 (Environmental)", "#27ae60"),
    "mono_Z3": ("Mono-Z3 (Coverage)", "#8e44ad"),
}
COMPROMISE_LABEL = ("Compromise (sol. 21)", "#2c5f8a")

LEVEL_ORDER = ["low", "base", "high"]

SCENARIOS = [
    ("price", "Price volatility",
     [("low", "price_0p70", "0.7x"), ("base", "base", "1.0x (base)"), ("high", "price_1p30", "1.3x")]),
    ("catchment", "Catchment sensitivity",
     [("low", "catchment_0p10", "0.10 km"), ("base", "base", "0.25 km (base)"), ("high", "catchment_0p50", "0.50 km")]),
    ("muncredit", "Municipal credit intensity",
     [("low", "muncredit_0p00", "0.00"), ("base", "base", "0.25 (base)"), ("high", "muncredit_0p50", "0.50")]),
    ("facility", "Facility cost/capacity",
     [("low", "facility_0p70", "0.7x"), ("base", "base", "1.0x (base)"), ("high", "facility_1p30", "1.3x")]),
]

plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 300,
})


def true_monthly_generation_kg():
    """Generacion REAL de toda la region (13 458 zonas), por mes, fija e
    independiente de cualquier escenario/nivel de sensibilidad."""
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


def series_for_all_levels(root_dir, run_ids, sol_number_fn, generation):
    """Return {level_key: pct_series} for whichever levels have data available
    (skips levels with no detalle_csv, e.g. an infeasible compromise run)."""
    out = {}
    for level_key, run_id in run_ids.items():
        detalle_dir = f"{root_dir}/{run_id}/detalle_csv"
        try:
            sol_number = sol_number_fn(detalle_dir)
            out[level_key] = coverage_pct_for_run(detalle_dir, sol_number, generation)
        except FileNotFoundError:
            continue
    return out


def main():
    generation = true_monthly_generation_kg()

    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), sharey=True)

    for ax, (sid, label, levels) in zip(axes, SCENARIOS):
        run_ids = {lvl: run for lvl, run, _desc in levels}

        all_solutions = list(MONO_SOL_LABELS.items()) + [(None, COMPROMISE_LABEL)]
        for estado, (sol_label, color) in all_solutions:
            if estado is None:
                by_level = series_for_all_levels(COMP_ROOT, run_ids, lambda d: 1, generation)
            else:
                by_level = series_for_all_levels(
                    MONO_ROOT, run_ids,
                    lambda d, estado=estado: mono_sol_map(d)[estado], generation)

            if "base" in by_level:
                base = by_level["base"]
                ax.plot(base.index, base.values, color=color, linewidth=2.0, zorder=3)

            if len(by_level) >= 2:
                stacked = pd.concat(by_level.values(), axis=1)
                lower, upper = stacked.min(axis=1), stacked.max(axis=1)
                ax.fill_between(lower.index, lower.values, upper.values,
                                 color=color, alpha=0.18, linewidth=0, zorder=1)

            if "low" not in by_level and estado is None:
                ax.text(0.5, 0.06, "Compromise: infeasible\nat lowest level",
                        transform=ax.transAxes, ha="center", va="bottom",
                        fontsize=7, color=COMPROMISE_LABEL[1])

        ax.set_title(label, fontsize=10)
        ax.set_xlabel("Month of year")
        ax.set_xticks(range(1, 13))

    axes[0].set_ylabel("% of regionwide\ngenerated waste collected")

    color_handles = [plt.Line2D([0], [0], color=c, lw=2, label=lbl)
                      for lbl, c in list(MONO_SOL_LABELS.values()) + [COMPROMISE_LABEL]]
    band_handle = [plt.Rectangle((0, 0), 1, 1, facecolor="gray", alpha=0.3,
                                   edgecolor="none", label="Range across low–high scenario levels")]

    fig.legend(handles=color_handles + band_handle, loc="lower center", ncol=5,
               frameon=False, bbox_to_anchor=(0.5, -0.06), fontsize=8.5)

    fig.suptitle("Monthly collection coverage vs. true regionwide waste generation\n"
                 "(solid line = base level; shaded band = range across the OFAT scenario)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PATH, bbox_inches="tight")
    print(f"Guardado: {OUT_PATH}")


if __name__ == "__main__":
    main()
