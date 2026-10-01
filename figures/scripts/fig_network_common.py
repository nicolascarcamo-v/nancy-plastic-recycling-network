"""Shared helpers to reconstruct the lost fig_network* scripts.

Data source: results/pareto_front/econ_results.json ("full" instance, 3111
residential zones I, 76 collection hubs J, 10 candidate facilities K).
Coordinates are UTM meters; converted to "local origin" km by subtracting the
bounding-box minimum of the residential zones (matches the axis ranges in the
original PDFs almost exactly, verified against known facility positions).

Backbone facilities = k in {1, 4, 7, 8}; any other opened k is an
"additional facility".
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # raiz del repositorio
import json

BACKBONE = {1, 4, 7, 8}

GOLD = "#d9a521"
MAROON = "#7a1f2b"
HUB_EDGE = "#1b3a6b"
EDGE_COLOR = "#5b7fa6"
BUILDING_COLOR = "0.55"


def load_geo(path=ROOT / "results/pareto_front/econ_results.json"):
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    geo = d["geo"]
    coords_i = {int(k): v for k, v in geo["coords_i"].items()}
    coords_j = {int(k): v for k, v in geo["coords_j"].items()}
    coords_k = {int(k): v for k, v in geo["coords_k"].items()}

    ox = min(v[0] for v in coords_i.values())
    oy = min(v[1] for v in coords_i.values())

    def to_km(x, y):
        return (x - ox) / 1000.0, (y - oy) / 1000.0

    zones = [to_km(v[0], v[1]) + (v[2],) for v in coords_i.values()]
    hubs = {j: to_km(v[0], v[1]) for j, v in coords_j.items()}
    facilities = {k: to_km(v[0], v[1]) + (v[2],) for k, v in coords_k.items()}
    return d, zones, hubs, facilities


def pareto_by_opened(pareto, opened):
    target = set(opened)
    for p in pareto:
        if set(p["opened"]) == target:
            return p
    raise KeyError(f"no pareto point with opened={opened}")


def sol_by_id(soluciones_rows, sol_id):
    for row in soluciones_rows:
        if int(row["sol"]) == sol_id:
            return row
    raise KeyError(f"sol {sol_id} not found")


def read_soluciones_csv(path=ROOT / "results/pareto_front/detalle_csv/soluciones.csv"):
    import csv
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def building_sizes(zones, s_min=3.0, s_max=260.0):
    w = [z[2] for z in zones]
    lo, hi = min(w) ** 0.5, max(w) ** 0.5
    return [s_min + (wi ** 0.5 - lo) / (hi - lo) * (s_max - s_min) for wi in w]


def draw_panel(ax, zones, hubs, facilities, opened, flows_jk, xlim, ylim,
                show_hub_marker=True):
    xi = [z[0] for z in zones]
    yi = [z[1] for z in zones]
    sizes = building_sizes(zones)
    ax.scatter(xi, yi, s=sizes, color=BUILDING_COLOR, alpha=0.35,
               linewidths=0, zorder=1)

    open_set = set(opened)
    for j, k, _flow in flows_jk:
        if k in open_set and j in hubs and k in facilities:
            x0, y0 = hubs[j]
            x1, y1 = facilities[k][0], facilities[k][1]
            ax.plot([x0, x1], [y0, y1], color=EDGE_COLOR,
                     linewidth=0.9, alpha=0.6, zorder=2)

    if show_hub_marker:
        xj = [v[0] for v in hubs.values()]
        yj = [v[1] for v in hubs.values()]
        ax.scatter(xj, yj, marker="s", s=22, facecolor="white",
                   edgecolor=HUB_EDGE, linewidths=0.8, zorder=3)

    for k, (x, y, name) in facilities.items():
        if k not in open_set:
            continue
        color = GOLD if k in BACKBONE else MAROON
        ax.scatter(x, y, marker="*", s=280, color=color,
                   edgecolor="k", linewidths=0.6, zorder=5)
        ax.annotate(f"k={k}", (x, y), textcoords="offset points",
                    xytext=(6, 4), fontsize=8, zorder=6)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.25, linewidth=0.6)


def shared_limits(zones, pad=0.35):
    xi = [z[0] for z in zones]
    yi = [z[1] for z in zones]
    return (min(xi) - pad, max(xi) + pad), (min(yi) - pad, max(yi) + pad)


def legend_handles():
    from matplotlib.lines import Line2D
    return [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=BUILDING_COLOR,
               alpha=0.6, markeredgecolor="none", markersize=10,
               label="Residential buildings"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="white",
               markeredgecolor=HUB_EDGE, markersize=8, label="Collection hubs"),
        Line2D([0], [0], marker="*", color="w", markerfacecolor=GOLD,
               markeredgecolor="k", markersize=14,
               label="Backbone facility ($k\\in\\{1,4,7,8\\}$)"),
        Line2D([0], [0], marker="*", color="w", markerfacecolor=MAROON,
               markeredgecolor="k", markersize=14, label="Additional facility"),
    ]
