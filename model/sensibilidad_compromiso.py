# -*- coding: utf-8 -*-
"""
sensibilidad_compromiso.py — análisis de sensibilidad de un factor (OFAT) sobre UNA
solución de compromiso concreta del frente de Pareto de la instancia Nancy.

Hermano de `sensibilidad_articulo.py`, que sensibiliza las tres soluciones
MONO-OBJETIVO (los anclas lexicográficos). Este script sensibiliza en cambio un
único punto del frente: la SOLUCIÓN DE COMPROMISO. Se mantienen como scripts
separados a propósito, para que quede auditoría de qué análisis produjo qué tabla.

Qué es la solución de compromiso y cómo se re-optimiza
------------------------------------------------------
El punto 21 del frente de `results/pareto_front` NO es el óptimo de ningún
objetivo aislado: es una celda de la grilla ε-constraint, es decir la solución de
un subproblema AUGMECON concreto (Mavrotas 2009, ec. (3) y (6)):

    max  Z4 + eps · ( s_Z2/r_Z2 + s_Z3/r_Z3 )
    s.a. Z2 + s_Z2 = eps_Z2        (eps_Z2 = -132,876.56)
         Z3 - s_Z3 = eps_Z3        (eps_Z3 =  60,035.88)
         s_k >= 0

Por eso "re-optimizar el punto" = **re-resolver ESE MISMO subproblema** bajo los
parámetros perturbados: mismo objetivo primario (Z4), mismas ε-constraints, mismos
RHS. Es un solve por configuración, no una tabla de pagos ni una grilla.

Los ε se mantienen FIJOS EN VALOR ABSOLUTO en todos los niveles. Es la elección
correcta para sensibilizar *una solución*: se fija la especificación (este nivel de
GEI y esta cobertura) y se perturba el parámetro, ceteris paribus. La alternativa
—recalcular la tabla de pagos en cada nivel y reposicionar el compromiso en el
mismo lugar RELATIVO del rango— responde otra pregunta y cuesta 9× más; se deja
disponible en `--eps-mode relativo` pero NO es el modo por defecto.

Consecuencia importante: con ε fijos, en algún nivel la especificación puede
volverse INFACTIBLE (p. ej. si al bajar el radio de captación ya no hay masa
suficiente para sostener Z3 >= 60,036). Eso NO es un fallo del script: es un
resultado, y de los más informativos ("por debajo de X km la solución de
compromiso deja de ser alcanzable"). Cuando ocurre, el script lanza un
diagnóstico que mide QUÉ ε es el inalcanzable y por cuánto — ver
`diagnostico_infactibles.csv`.

Punto de referencia (results/pareto_front, punto 21 del frente)
---------------------------------------------------------------
    n_open   7   ->  plantas [1, 4, 7, 8, 9, 10, 14]
    Z4 (VAN)      146,887.81 EUR
    Z2 (GEI)     -142,544.88 kgCO2eq
    Z3 (cobert.)   60,035.71
    eps_Z2       -132,876.56      eps_Z3   60,035.88
    status       timelimit,  gap 5.80%

⚠ El punto de referencia cerró con 5.80% de gap (topó el límite de 1000 s): es un
incumbente, no un óptimo probado. Los deltas contra él heredan esa incertidumbre.
La columna gap_pct reporta lo mismo en cada nivel perturbado.

La solución de referencia se usa además como WARM START en cada nivel.

Escenarios (tabla del manuscrito) — idénticos a sensibilidad_articulo.py
-------------------------------------------------------------------------
  Price volatility            PRICE_FACTOR             1.0    ×{0.7, 1.0, 1.3}
  Catchment sensitivity       D_BAR_IJ                 0.25   {0.10, 0.25, 0.50} km
  Municipal credit intensity  MUNICIPAL_CREDIT_SHARE   0.25   {0.00, 0.25, 0.50}
  Facility cost/capacity      F_fix y Q_prod (juntos)  1.0    ×{0.7, 1.0, 1.3}

9 configuraciones únicas (el nivel base se comparte entre los cuatro escenarios)
× 1 solve = 9 solves.

USO
---
  python sensibilidad_compromiso.py --dry-run            # plan + tamaños, sin resolver
  python sensibilidad_compromiso.py                      # los 4 escenarios, zona full
  python sensibilidad_compromiso.py --zone centro        # prueba rápida
  python sensibilidad_compromiso.py --only catchment --outdir results/sensitivity_compromise
  python sensibilidad_compromiso.py --point 21           # otro punto del frente
  python sensibilidad_compromiso.py --force              # re-resolver lo ya hecho

Reanuda por defecto: una configuración con `compromiso_results.json` ya escrito se
salta y se relee. `--force` la re-resuelve.

SALIDAS (bajo --outdir, por defecto results_compromiso_<zona>_<timestamp>/)
----------------------------------------------------------------------------
  runs/<run_id>/compromiso_results.json  solución completa + config + ε usados
  runs/<run_id>/detalle_csv/*.csv        plantas, hubs, flujos, zonas
  runs/<run_id>/run_config.json          configuración exacta de esa corrida
  runs/<run_id>/run.log                  stdout completo del solve
  sensibilidad_compromiso.csv    1 fila por (escenario, nivel) — TABLA PRINCIPAL
  cardinalidad_compromiso.csv    n_plantas y qué plantas entran/salen vs el punto 21
  diagnostico_infactibles.csv    niveles donde la especificación ε no es alcanzable
  plantas_detalle.csv            1 fila por (escenario, nivel, planta abierta)
  runs_manifest.json             config + estado + tiempos + ε de cada corrida
  comparacion_*.png              cardinalidad y VAN por nivel, con marcas de infactible

Las tablas se reescriben DESPUÉS DE CADA CORRIDA: si se interrumpe, lo resuelto
queda en disco y es legible.
"""

import os
import sys
import csv
import json
import time
import argparse
import datetime
import traceback

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")   # los SystemExit salen por aquí
except Exception:
    pass

import pyomo.environ as pyo
import sensibilidad as sens


# ===========================================================================
# CONFIGURACIÓN
# ===========================================================================

# Corrida y punto de referencia: de aquí sale la solución de compromiso que se
# sensibiliza (y que sirve de warm start y de término de comparación).
REFERENCE_DIR = os.path.join(sens.REPO_ROOT, "results", "pareto_front")
REFERENCE_POINT = 21          # índice 1-based dentro de payload["pareto"]

# Subproblema AUGMECON que define el punto: primario + objetivos con ε-constraint.
# Debe coincidir con la configuración que generó el frente de referencia
# (ECON_PRIMARY="Z4", ECON_CONSTRAINED_ALL=["Z2","Z3"]).
PRIMARY = "Z4"
CONSTRAINED = ["Z2", "Z3"]

# Un solo solve por configuración: se puede dar más presupuesto que en un barrido
# de 81 solves. El punto de referencia se resolvió con TL=1000 y quedó en 5.80% de
# gap, así que 1000 s es el mínimo razonable para que los niveles sean comparables
# con él y entre sí.
TIME_LIMIT_DEFAULT = 1000
MIP_GAP_DEFAULT = 0.005
ZONE_DEFAULT = "full"

# Caso base = escenario BASE de sensibilidad.py + el factor de plantas que este
# script añade. Coincide con los params de la corrida de referencia.
BASE_CFG = dict(
    D_BAR_IJ=0.25,
    D_BAR_JK=1.5,
    D_BAR_KK_KM=0.6,
    D_BAR_CENTER_K=None,
    N_CAP=None,
    BETA_MIN=0.0,
    WASTE_FACTOR=1.0,
    LCA_CREDIT=True,
    EPS_TR_SCALE=1.0,
    PRICE_FACTOR=1.0,
    MUNICIPAL_CREDIT_SHARE=0.25,
    FACILITY_FACTOR=1.0,      # escala F_k y Q̄^prod_k CONJUNTAMENTE
)

#   (id, etiqueta del manuscrito, clave de BASE_CFG que se barre, niveles)
SWEEPS = [
    ("price",     "Price volatility",           "PRICE_FACTOR",           [0.7, 1.0, 1.3]),
    ("catchment", "Catchment sensitivity",      "D_BAR_IJ",               [0.10, 0.25, 0.50]),
    ("muncredit", "Municipal credit intensity", "MUNICIPAL_CREDIT_SHARE", [0.00, 0.25, 0.50]),
    ("facility",  "Facility cost/capacity",     "FACILITY_FACTOR",        [0.7, 1.0, 1.3]),
    ("beta",      "Minimum collection floor",   "BETA_MIN",               [0.0, 0.05, 0.10, 0.20]),
]


# ===========================================================================
# UTILIDADES
# ===========================================================================
class _Tee:
    """Duplica stdout a un archivo, para que cada corrida deje su propio log."""

    def __init__(self, stream, path):
        self.stream = stream
        self.fh = open(path, "w", encoding="utf-8")

    def write(self, s):
        self.stream.write(s)
        self.fh.write(s)

    def flush(self):
        self.stream.flush()
        self.fh.flush()

    def close(self):
        try:
            self.fh.close()
        except Exception:
            pass


def _fmt_level(v):
    """Nombre de carpeta estable para un nivel numérico: 0.7 -> '0p70'."""
    return f"{v:.2f}".replace(".", "p").replace("-", "m")


def check_instance_cache(zone, cfg_list):
    """El caché sólo guarda arcos hasta arc_ij_max / arc_jk_max. Un D_BAR_IJ por
    encima de ese techo NO amplía el grafo: devolvería el mismo grafo recortado y
    el barrido de catchment saldría plano por artefacto del caché, no por el
    modelo."""
    man_path = os.path.join(sens.INSTANCES_DIR, "manifest.json")
    if not os.path.exists(man_path):
        print(f"  [warning] Sin {man_path}: no puedo verificar el techo de arcos del caché.")
        return
    with open(man_path, encoding="utf-8") as fh:
        zones = json.load(fh).get("zones", [])
    z = next((x for x in zones if x.get("name") == zone), None)
    if z is None:
        print(f"  [warning] La zona '{zone}' no está en el manifest; no verifico arcos.")
        return

    ij_max = float(z.get("arc_ij_max", 0.0))
    jk_max = float(z.get("arc_jk_max", 0.0))
    need_ij = max(c["D_BAR_IJ"] for c in cfg_list)
    need_jk = max(c["D_BAR_JK"] for c in cfg_list)

    if need_ij > ij_max + 1e-9:
        raise SystemExit(
            f"\nD_BAR_IJ máximo del barrido = {need_ij} km, pero el caché de '{zone}' sólo "
            f"guarda arcos I-J hasta {ij_max} km.\nRegenera el caché (build_instances.py, "
            f"ajustando arc_ij_max): si no, el nivel alto de 'catchment' sería idéntico "
            f"al techo del caché y el escenario saldría plano por artefacto.")
    if need_jk > jk_max + 1e-9:
        raise SystemExit(
            f"\nD_BAR_JK = {need_jk} km supera el techo del caché ({jk_max} km). "
            f"Regenera con build_instances.py.")

    if abs(need_ij - ij_max) < 1e-9:
        print(f"  [info] El nivel alto de catchment ({need_ij} km) coincide EXACTAMENTE con "
              f"el techo del caché (arc_ij_max={ij_max}). Válido, pero no se puede subir "
              f"más sin re-correr build_instances.py.")


# ===========================================================================
# REFERENCIA: LA SOLUCIÓN DE COMPROMISO
# ===========================================================================
def load_reference_point(refdir, index):
    """Lee el punto `index` (1-based) del frente de la corrida de referencia y la
    tabla de pagos que definió las escalas r_k del término aumentado.

    Devuelve (punto, eps, ranges, params). `eps` son los RHS de las ε-constraints
    que generaron el punto; `ranges` es el rango [min,max] por objetivo de la
    tabla de pagos, necesario para reconstruir EL MISMO objetivo aumentado."""
    path = os.path.join(refdir, "econ_results.json")
    if not os.path.isabs(path):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
    if not os.path.exists(path):
        raise SystemExit(f"\nNo existe la corrida de referencia: {path}\n"
                         f"Ajusta REFERENCE_DIR o pasa --reference <carpeta>.")
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)

    front = data.get("pareto", [])
    if not front:
        raise SystemExit(f"\n{path} no tiene frente de Pareto ('pareto' vacío).")
    if not (1 <= index <= len(front)):
        raise SystemExit(f"\nEl punto {index} está fuera de rango: el frente de {path} "
                         f"tiene {len(front)} puntos (1..{len(front)}).")

    pt = front[index - 1]

    eps = {}
    for k in CONSTRAINED:
        key = f"eps_{k}"
        if key not in pt:
            raise SystemExit(
                f"\nEl punto {index} no tiene '{key}': no proviene de una celda de la "
                f"grilla ε-constraint (¿es un ancla del payoff?). Este script sensibiliza "
                f"una solución de COMPROMISO, que por definición sale de la grilla. "
                f"Para sensibilizar los anclas mono-objetivo usa sensibilidad_articulo.py.")
        eps[k] = float(pt[key])

    ranges = sens._rng_from_payoff(data.get("payoff", {}), CONSTRAINED)
    missing = [k for k in CONSTRAINED if k not in ranges]
    if missing:
        raise SystemExit(f"\nLa tabla de pagos de la referencia no tiene fila(s) para "
                         f"{missing}: no puedo reconstruir las escalas r_k del término "
                         f"aumentado sin ellas.")

    print(f"  referencia : {path}")
    print(f"  punto      : {index} de {len(front)} (solución de compromiso)")
    gp = pt.get("gap_pct")
    gtxt = f"{gp:.2f}%" if isinstance(gp, (int, float)) else str(gp)
    warn = "   <-- incumbente, NO óptimo probado" if isinstance(gp, (int, float)) and gp > 1 else ""
    print(f"    n_open={pt.get('n_open')}  plantas={pt.get('opened')}")
    print(f"    Z4={pt.get('Z4', 0):,.1f}   Z2={pt.get('Z2', 0):,.1f}   Z3={pt.get('Z3', 0):,.1f}")
    print(f"    status={pt.get('status')}  gap={gtxt}{warn}")
    print(f"    ε fijos del subproblema: " + "   ".join(f"{k}={v:,.2f}" for k, v in eps.items()))
    print(f"    rangos r_k del payoff  : "
          + "   ".join(f"{k}=[{ranges[k][0]:,.1f}, {ranges[k][1]:,.1f}]" for k in CONSTRAINED))
    return pt, eps, ranges, data.get("params", {})


def seed_from_reference(m, D, ref_pt):
    """Siembra y[k,t] con la solución de compromiso, como MIP start.

    Sólo se siembran plantas que EXISTEN en la instancia perturbada: al mover
    d̄^IJ el conjunto K cambia, así que una planta del punto de referencia puede
    no estar presente. Es un arranque PARCIAL (sólo binarias); el solver completa
    el resto, y si el backend no lo acepta `_safe_solve` reintenta sin warmstart."""
    if not ref_pt:
        return 0
    ty = ref_pt.get("temporal_y", {}) or {}
    K, TT = set(D["K"]), set(D["TT"])
    seeded = 0
    for k in D["K"]:                                   # arranca todo cerrado
        for t in D["TT"]:
            m.y[k, t].value = 0
    for k_str, serie in ty.items():
        k = int(k_str)
        if k not in K:
            continue
        for t_str, v in serie.items():
            t = int(t_str)
            if t not in TT:
                continue
            m.y[k, t].value = int(v)
            seeded += 1
    return seeded


# ===========================================================================
# APLICACIÓN DE LA CONFIGURACIÓN
# ===========================================================================
def apply_config(cfg, run_id, time_limit, mip_gap):
    """Escribe la configuración en los globals de sensibilidad.py. Funciona
    porque `load_instance`, `_build_final_model_dict` y `build_model` los leen EN
    TIEMPO DE LLAMADA, no al importar."""
    sens.D_BAR_IJ = cfg["D_BAR_IJ"]
    sens.D_BAR_JK = cfg["D_BAR_JK"]
    sens.D_BAR_KK_KM = cfg["D_BAR_KK_KM"]
    sens.D_BAR_CENTER_K = cfg["D_BAR_CENTER_K"]
    sens.N_CAP = cfg["N_CAP"]
    sens.BETA_MIN = cfg["BETA_MIN"]
    sens.WASTE_FACTOR = cfg["WASTE_FACTOR"]
    sens.LCA_CREDIT = cfg["LCA_CREDIT"]
    sens.EPS_TR_SCALE = cfg["EPS_TR_SCALE"]
    sens.PRICE_FACTOR = cfg["PRICE_FACTOR"]
    sens.MUNICIPAL_CREDIT_SHARE = cfg["MUNICIPAL_CREDIT_SHARE"]

    sens.SCENARIO = run_id
    sens.SOLVER_TIME_LIMIT = time_limit      # make_solver() lo lee: fijarlo ANTES
    sens.SOLVER_MIP_GAP = mip_gap
    sens.ECON_PRIMARY = PRIMARY
    sens.ECON_CONSTRAINED_ALL = list(CONSTRAINED)


def scale_facilities(D, factor):
    """Escenario 4: escala F_k y Q̄^prod_k CONJUNTAMENTE, como pide el manuscrito.

    No hay global para esto: F_fix y Q_prod se leen directo del caché. Se escalan
    sobre el dict D ya construido, que es lo que consume `build_model`:
        · D["Q_prod"][k] -> restricción intake_cap (Σx ≤ Q_prod·y)
        · D["F_fix"][k]  -> costo fijo en π_k
    Son sus dos únicos consumidores, así que equivale a escalar la fuente.

    M_pi_global (big-M de Z1) se recalcula porque depende de max(Q_prod): si se
    escala hacia arriba sin tocarlo, el big-M queda POR DEBAJO del π_k alcanzable
    y phi_constraint recortaría soluciones factibles."""
    if abs(factor - 1.0) < 1e-12:
        return D

    D["Q_prod"] = {k: v * factor for k, v in D["Q_prod"].items()}
    D["F_fix"] = {k: v * factor for k, v in D["F_fix"].items()}

    disc_sum = sum(D["disc"].values())
    max_margin = max(
        (D["price"][mm] * D["eta"][mm] + sens.C_MUN * D["municipal_credit_share"])
        for mm in D["M"]
    ) if D["M"] else 0.0
    max_Qprod = max(D["Q_prod"].values()) if D["Q_prod"] else 0.0
    big_m = 1.05 * disc_sum * max_Qprod * max_margin
    D["M_pi_global"] = big_m if big_m > 0 else 1e6

    print(f"  [info] Plantas escaladas ×{factor}: "
          f"F_k ∈ [{min(D['F_fix'].values()):,.0f}, {max(D['F_fix'].values()):,.0f}] EUR, "
          f"Q̄^prod ∈ [{min(D['Q_prod'].values()):,.0f}, {max(D['Q_prod'].values()):,.0f}] kg, "
          f"big-M = {D['M_pi_global']:,.0f}")
    return D


# ===========================================================================
# RE-OPTIMIZACIÓN DE LA SOLUCIÓN DE COMPROMISO
# ===========================================================================
def _loose_eps(m, R, k, ranges):
    """Deja la ε-constraint de k no vinculante (en su nadir)."""
    lo, hi = ranges[k]
    sens._set_eps(m, R, k, lo if sens.OBJ_SENSE[k] == pyo.maximize else hi)


def diagnose_infeasibility(m, solver, D, eps, R, ranges):
    """Cuando la especificación ε del compromiso resulta infactible, mide QUÉ
    restricción es la inalcanzable y por cuánto.

    Para cada objetivo con ε: se afloja el OTRO y se optimiza éste en solitario,
    obteniendo su mejor valor alcanzable bajo los parámetros perturbados. Si ese
    óptimo no llega al ε requerido, esa es la restricción que rompe el punto.
    Son 2 solves y sólo se pagan cuando ya hubo una infactibilidad."""
    print("    diagnóstico de infactibilidad (qué ε no se alcanza):")
    diag = {}
    for k in CONSTRAINED:
        for other in CONSTRAINED:                      # afloja todos...
            _loose_eps(m, R, other, ranges)
        sens.set_objective(m, sens._obj_expr(m, k), sens.OBJ_SENSE[k])
        cat, _, meta = sens._safe_solve(solver, m, tee=sens.SOLVER_TEE, warmstart=False)
        if cat in ("infeasible", "nosolution", "error"):
            diag[k] = dict(alcanzable=None, requerido=eps[k], cumple=None, estado=cat)
            print(f"      {k}: {cat} — no se pudo medir el óptimo aislado")
            continue
        best = float(pyo.value(sens._obj_expr(m, k)))
        if sens.OBJ_SENSE[k] == pyo.maximize:
            ok = best >= eps[k] - 1e-6
            holgura = best - eps[k]
        else:
            ok = best <= eps[k] + 1e-6
            holgura = eps[k] - best
        diag[k] = dict(alcanzable=best, requerido=eps[k], cumple=bool(ok),
                       holgura=holgura, estado=cat, gap_pct=meta["gap_pct"])
        sent = "≥" if sens.OBJ_SENSE[k] == pyo.maximize else "≤"
        print(f"      {k}: mejor alcanzable={best:>14,.1f}   requerido {sent}{eps[k]:>14,.1f}   "
              f"{'CUMPLE' if ok else 'NO ALCANZA'}  (holgura={holgura:,.1f})")
    return diag


def solve_compromise(m, solver, D, ref_pt, eps, ranges, augw):
    """Re-resuelve el subproblema AUGMECON que define la solución de compromiso,
    bajo la configuración de parámetros vigente.

    Reconstruye exactamente el subproblema del punto de referencia:
      · ε-constraints como IGUALDADES con holgura (Mavrotas 2009, ec. (3)) vía
        `_ensure_augmecon`;
      · objetivo aumentado max Z4 + eps·Σ s_k/r_k (ec. (6)) vía
        `_setup_augmecon_objective`;
      · r_k tomados de la tabla de pagos DE LA REFERENCIA, no recalculados: son
        parte de la definición del subproblema que generó el punto, y recalcularlos
        cambiaría el objetivo entre niveles y rompería la comparabilidad;
      · RHS ε fijos en el valor absoluto del punto de referencia.

    Devuelve (punto | None, meta, diagnóstico | None)."""
    sens._ensure_augmecon(m, CONSTRAINED)
    R = sens._augmecon_ranges(ranges, CONSTRAINED)
    sens._setup_augmecon_objective(m, PRIMARY, CONSTRAINED, R, augw=augw)
    for k in CONSTRAINED:
        sens._set_eps(m, R, k, eps[k])

    n_seed = seed_from_reference(m, D, ref_pt)
    rel = {k: ("≤" if sens.OBJ_SENSE[k] == pyo.minimize else "≥") for k in CONSTRAINED}
    print(f"    subproblema: max {PRIMARY} + {augw}·Σ s_k/r_k   s.a.  "
          + "   ".join(f"{k}{rel[k]}{eps[k]:,.1f}" for k in CONSTRAINED))
    if n_seed:
        print(f"    warm start: {n_seed} binarias desde la solución de referencia")

    cat, _, meta = sens._safe_solve(solver, m, tee=sens.SOLVER_TEE, warmstart=True)
    gp = meta["gap_pct"]
    gtxt = f"{gp:.2f}%" if isinstance(gp, (int, float)) else str(gp)
    print(f"    -> {cat}   gap={gtxt}   t={meta['wall_time_s']:.0f}s")

    if cat == "infeasible":
        # La especificación del compromiso no es alcanzable a este nivel. Es un
        # resultado, no un error: se mide cuál de los dos ε la rompe.
        diag = diagnose_infeasibility(m, solver, D, eps, R, ranges)
        return None, meta, diag

    if cat in ("nosolution", "error"):
        print(f"    (sin incumbente en {sens.SOLVER_TIME_LIMIT}s; no es infactibilidad "
              f"probada — subir --time-limit)")
        return None, meta, None

    pt = sens._record_point(m, status="compromiso", termination=cat,
                            gap_pct=meta["gap_pct"], bound=meta["bound"],
                            wall_time_s=meta["wall_time_s"],
                            **{f"eps_{k}": float(v) for k, v in eps.items()})
    return pt, meta, None


# ===========================================================================
# PLAN DE CORRIDAS
# ===========================================================================
def build_plan():
    """Devuelve (runs, plan). El nivel que coincide con el base apunta a la única
    corrida `base`, compartida por los cuatro escenarios."""
    runs = {"base": dict(BASE_CFG)}
    plan = [("base", "Base", "—", None, "base", True)]

    for sid, label, param, levels in SWEEPS:
        for lv in levels:
            cfg = dict(BASE_CFG)
            cfg[param] = lv
            is_base = abs(lv - BASE_CFG[param]) < 1e-12
            run_id = "base" if is_base else f"{sid}_{_fmt_level(lv)}"
            if run_id not in runs:
                runs[run_id] = cfg
            plan.append((sid, label, param, lv, run_id, is_base))
    return runs, plan


# ===========================================================================
# UNA CORRIDA
# ===========================================================================
def solve_run(run_id, cfg, zone, rundir, ref_pt, eps, ranges,
              time_limit, mip_gap, augw, point_index):
    """Resuelve la solución de compromiso bajo una configuración. Nunca lanza:
    los errores se capturan para que una corrida caída no aborte el barrido."""
    os.makedirs(rundir, exist_ok=True)

    cfg_payload = dict(
        run_id=run_id, zone=zone, config=cfg,
        reference=REFERENCE_DIR, reference_point=point_index,
        subproblem=dict(primary=PRIMARY, constrained=CONSTRAINED, eps=eps,
                        payoff_ranges={k: list(v) for k, v in ranges.items()}, augw=augw),
        solver=dict(time_limit_s=time_limit, mip_gap=mip_gap),
        horizon=sens.HORIZON, alpha=sens.ALPHA,
        yield_mode=sens.YIELD_MODE, use_disposal=sens.USE_DISPOSAL,
        lca_credit_m=sens.LCA_CREDIT_M, c_mun=sens.C_MUN,
        augmecon=dict(bypass=sens.AUGMECON2_BYPASS,
                      eps_feas_tol_frac=sens.EPS_FEAS_TOL_FRAC,
                      tiebreak_open_w=sens.TIEBREAK_OPEN_W),
        started_at=datetime.datetime.now().isoformat(timespec="seconds"),
    )
    with open(os.path.join(rundir, "run_config.json"), "w", encoding="utf-8") as fh:
        json.dump(cfg_payload, fh, indent=1, ensure_ascii=False)

    tee = _Tee(sys.stdout, os.path.join(rundir, "run.log"))
    real_stdout = sys.stdout
    sys.stdout = tee
    t0 = time.time()
    status = dict(run_id=run_id, status="ok", error=None, factible=True)
    punto, diag = None, None

    try:
        apply_config(cfg, run_id, time_limit, mip_gap)

        print(f"\nCargando zona '{zone}' desde caché ({sens.INSTANCES_DIR}) ...")
        D = sens.load_instance(zone)
        if D is None:
            raise RuntimeError(f"load_instance devolvió None para la zona '{zone}'.")
        scale_facilities(D, cfg["FACILITY_FACTOR"])

        print(f"  |I|={len(D['I'])}  |J|={len(D['J'])}  |K|={len(D['K'])}  "
              f"|M|={len(D['M'])}  |T|={len(D['TT'])}  "
              f"|A_IJ|={len(D['A_IJ'])}  |A_JK|={len(D['A_JK'])}")
        status.update(n_I=len(D["I"]), n_J=len(D["J"]), n_K=len(D["K"]),
                      n_arc_ij=len(D["A_IJ"]), n_arc_jk=len(D["A_JK"]))

        tb = time.time()
        m = sens.build_model(D)
        nv = sum(1 for _ in m.component_data_objects(pyo.Var, active=True))
        nc = sum(1 for _ in m.component_data_objects(pyo.Constraint, active=True))
        print(f"  modelo en {time.time()-tb:.1f}s: {nv:,} vars "
              f"({len(D['K'])*len(D['TT']):,} binarias), {nc:,} restricciones")
        status.update(n_vars=nv, n_cons=nc)

        solver = sens.make_solver()

        print(f"\n  == solución de compromiso (punto {point_index} de la referencia) ==")
        punto, meta, diag = solve_compromise(m, solver, D, ref_pt, eps, ranges, augw)

        if punto is None:
            status["factible"] = False
            status["status"] = "infactible" if diag is not None else "sin_solucion"
            print(f"    -> la especificación del compromiso NO se sostiene en este nivel")
        else:
            print(f"    -> n_open={punto['n_open']}  plantas={punto['opened']}  "
                  f"Z4={punto.get('Z4', 0):,.0f}  Z2={punto.get('Z2', 0):,.0f}  "
                  f"Z3={punto.get('Z3', 0):,.1f}")

        save_run(punto, diag, D, cfg, run_id, zone, rundir, eps, ranges, point_index)

    except Exception as exc:
        status["status"] = "error"
        status["error"] = f"{type(exc).__name__}: {exc}"
        print(f"\n  [ERROR] La corrida '{run_id}' falló: {status['error']}")
        traceback.print_exc(file=sys.stdout)
    finally:
        status["wall_time_s"] = round(time.time() - t0, 1)
        status["finished_at"] = datetime.datetime.now().isoformat(timespec="seconds")
        sys.stdout = real_stdout
        tee.close()

    cfg_payload["result"] = status
    with open(os.path.join(rundir, "run_config.json"), "w", encoding="utf-8") as fh:
        json.dump(cfg_payload, fh, indent=1, ensure_ascii=False)
    return status


def save_run(punto, diag, D, cfg, run_id, zone, rundir, eps, ranges, point_index):
    """Guarda la solución de una configuración: JSON completo + los CSV detallados
    de sensibilidad.py (operan sobre puntos de `_record_point`, así que sirven
    igual para una solución única)."""
    payload = dict(
        run_id=run_id, zone=zone, config=cfg,
        reference=REFERENCE_DIR, reference_point=point_index,
        subproblem=dict(primary=PRIMARY, constrained=CONSTRAINED, eps=eps,
                        payoff_ranges={k: list(v) for k, v in ranges.items()}),
        factible=punto is not None,
        solucion=sens._strip_heavy(punto, keep_full=True) if punto else None,
        diagnostico_infactibilidad=diag,
        geo=sens._build_geo(D),
    )
    path = os.path.join(rundir, "compromiso_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False)
    print(f"\n  resultados → {path}")

    if punto is not None:
        try:
            sens.save_detailed_csv([punto], D, rundir)
        except Exception as exc:
            print(f"  [warning] save_detailed_csv falló: {exc}")


def load_run(rundir):
    """Relee una corrida ya resuelta (modo reanudar)."""
    path = os.path.join(rundir, "compromiso_results.json")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


# ===========================================================================
# TABLAS AGREGADAS
# ===========================================================================
def _set_delta(ref_opened, new_opened):
    """Qué plantas entran y salen respecto de la solución de compromiso."""
    r = set(int(x) for x in (ref_opened or []))
    n = set(int(x) for x in (new_opened or []))
    entran, salen = sorted(n - r), sorted(r - n)
    jac = (len(r & n) / len(r | n)) if (r | n) else 1.0
    return entran, salen, round(jac, 4)


def write_aggregates(outdir, plan, results, runs, ref_pt, eps, point_index):
    """Reescribe las tablas agregadas + el manifest. Se llama tras CADA corrida."""

    def _w(name, header, rows):
        with open(os.path.join(outdir, name), "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(header)
            wr.writerows(rows)

    ref_open = ref_pt.get("opened", [])
    ref_n = ref_pt.get("n_open")

    # --- 1. TABLA PRINCIPAL ---------------------------------------------------
    hdr = ["escenario", "etiqueta", "parametro", "nivel", "es_base", "run_id",
           "estado", "factible",
           "Z1_VAN_peor_EUR", "Z2_GHG_kgCO2e", "Z3_cobertura", "Z4_NPV_EUR", "Z5_output_kg",
           "cobertura_pct", "colectado_total_kg",
           "n_plantas", "plantas",
           "ref_n_plantas", "ref_plantas", "plantas_entran", "plantas_salen", "jaccard_vs_ref",
           "delta_Z4_vs_ref", "delta_Z4_pct", "delta_Z2_vs_ref", "delta_Z3_vs_ref",
           "eps_Z2", "eps_Z3", "gap_pct", "bound", "wall_time_s", "terminacion", "error"]
    rows = []
    for sid, label, param, lv, run_id, is_base in plan:
        r = results.get(run_id, {})
        st = r.get("status", {})
        p = r.get("solucion")
        base_cols = [sid, label, param, ("" if lv is None else lv), int(is_base), run_id,
                     st.get("status", "pendiente"), int(bool(r.get("factible", False)))]
        if not p:
            rows.append(base_cols + [""] * 23 + [st.get("error", "") or ""])
            continue
        entran, salen, jac = _set_delta(ref_open, p.get("opened"))
        z4, z2, z3 = p.get("Z4", 0.0), p.get("Z2", 0.0), p.get("Z3", 0.0)
        rz4, rz2, rz3 = ref_pt.get("Z4", 0.0), ref_pt.get("Z2", 0.0), ref_pt.get("Z3", 0.0)
        rows.append(base_cols + [
            p.get("Z1", ""), z2, z3, z4, p.get("Z5", ""),
            p.get("coverage_pct", ""), p.get("collected_total", ""),
            p.get("n_open", ""), " ".join(map(str, p.get("opened", []))),
            ref_n, " ".join(map(str, ref_open)),
            " ".join(map(str, entran)), " ".join(map(str, salen)), jac,
            z4 - rz4, (100.0 * (z4 - rz4) / abs(rz4)) if rz4 else "",
            z2 - rz2, z3 - rz3,
            p.get("eps_Z2", ""), p.get("eps_Z3", ""),
            p.get("gap_pct", "NA"), p.get("bound", ""),
            p.get("wall_time_s", ""), p.get("termination", ""), "",
        ])
    _w("sensibilidad_compromiso.csv", hdr, rows)

    # --- 2. CARDINALIDAD ------------------------------------------------------
    hdr = ["escenario", "parametro", "nivel", "factible",
           "n_plantas", "ref_n_plantas", "delta_n_plantas",
           "plantas", "plantas_entran", "plantas_salen", "jaccard_vs_ref",
           "misma_cardinalidad", "misma_identidad", "gap_pct"]
    rows = []
    for sid, label, param, lv, run_id, is_base in plan:
        r = results.get(run_id, {})
        p = r.get("solucion")
        if not p:
            rows.append([sid, param, ("" if lv is None else lv), 0,
                         "", ref_n, "", "", "", "", "", "", "", ""])
            continue
        entran, salen, jac = _set_delta(ref_open, p.get("opened"))
        n = p.get("n_open")
        rows.append([sid, param, ("" if lv is None else lv), 1,
                     n, ref_n,
                     (n - ref_n) if isinstance(n, int) and isinstance(ref_n, int) else "",
                     " ".join(map(str, p.get("opened", []))),
                     " ".join(map(str, entran)), " ".join(map(str, salen)), jac,
                     int(n == ref_n) if isinstance(n, int) and isinstance(ref_n, int) else "",
                     int(not entran and not salen),
                     p.get("gap_pct", "NA")])
    _w("cardinalidad_compromiso.csv", hdr, rows)

    # --- 3. DIAGNÓSTICO DE INFACTIBILIDAD -------------------------------------
    hdr = ["escenario", "parametro", "nivel", "run_id", "objetivo", "sentido",
           "requerido_eps", "mejor_alcanzable", "cumple", "holgura", "estado_solve"]
    rows = []
    for sid, label, param, lv, run_id, is_base in plan:
        diag = results.get(run_id, {}).get("diagnostico_infactibilidad")
        if not diag:
            continue
        for k, dd in diag.items():
            rows.append([sid, param, ("" if lv is None else lv), run_id, k,
                         "max" if sens.OBJ_SENSE[k] == pyo.maximize else "min",
                         dd.get("requerido", ""), dd.get("alcanzable", ""),
                         ("" if dd.get("cumple") is None else int(dd["cumple"])),
                         dd.get("holgura", ""), dd.get("estado", "")])
    _w("diagnostico_infactibles.csv", hdr, rows)

    # --- 4. DETALLE POR PLANTA ------------------------------------------------
    hdr = ["escenario", "nivel", "run_id", "k", "nombre", "x", "y",
           "apertura_t", "VAN_planta_EUR", "en_referencia"]
    rows = []
    refset = set(int(x) for x in ref_open)
    for sid, label, param, lv, run_id, is_base in plan:
        r = results.get(run_id, {})
        p = r.get("solucion")
        if not p:
            continue
        coords = r.get("geo", {}).get("coords_k", {})
        pi = p.get("pi_individual", {}) or {}
        for k in p.get("opened", []):
            ck = coords.get(str(k), [None, None, ""])
            rows.append([sid, ("" if lv is None else lv), run_id, k,
                         (ck[2] if len(ck) > 2 else ""), ck[0], ck[1],
                         p.get("open_period", {}).get(str(k),
                               p.get("open_period", {}).get(k, "")),
                         pi.get(str(k), pi.get(k, "")),
                         int(int(k) in refset)])
    _w("plantas_detalle.csv", hdr, rows)

    # --- 5. manifest -----------------------------------------------------------
    manifest = dict(
        generated_at=datetime.datetime.now().isoformat(timespec="seconds"),
        analysis="one-factor (OFAT) sensitivity of a single compromise solution; the "
                 "AUGMECON subproblem that generated it is re-solved at each level "
                 "with the epsilon RHS held fixed in absolute value",
        reference=REFERENCE_DIR, reference_point=point_index,
        reference_solution=dict(n_open=ref_n, opened=ref_open,
                                Z4=ref_pt.get("Z4"), Z2=ref_pt.get("Z2"),
                                Z3=ref_pt.get("Z3"), gap_pct=ref_pt.get("gap_pct"),
                                status=ref_pt.get("status")),
        subproblem=dict(primary=PRIMARY, constrained=CONSTRAINED, eps=eps),
        base_config=BASE_CFG,
        sweeps=[dict(id=s, label=l, param=p, levels=v) for s, l, p, v in SWEEPS],
        runs={rid: dict(config=cfg,
                        status=results.get(rid, {}).get("status", {}),
                        factible=results.get(rid, {}).get("factible"))
              for rid, cfg in runs.items()},
        plan=[dict(escenario=s, etiqueta=l, parametro=p, nivel=lv, run_id=r, es_base=b)
              for s, l, p, lv, r, b in plan],
    )
    with open(os.path.join(outdir, "runs_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1, ensure_ascii=False)


# ===========================================================================
# GRÁFICOS
# ===========================================================================
def make_comparison_plots(outdir, plan, results, ref_pt):
    """Por escenario: cardinalidad y VAN por nivel, con los niveles infactibles
    marcados explícitamente (no se omiten: que la especificación no se sostenga
    es parte del resultado)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print("  matplotlib no disponible; omito los gráficos comparativos.")
        return

    ref_n = ref_pt.get("n_open")
    ref_z4 = ref_pt.get("Z4")

    for sid, label, param, levels in SWEEPS:
        entries = [(lv, results.get(rid, {})) for s, l, p, lv, rid, b in plan if s == sid]
        if not entries:
            continue

        ok = [(lv, r["solucion"]) for lv, r in entries if r.get("solucion")]
        bad = [lv for lv, r in entries
               if not r.get("solucion") and r.get("status", {}).get("status") == "infactible"]
        if not ok and not bad:
            continue

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

        if ok:
            xs = [lv for lv, _ in ok]
            axes[0].plot(xs, [p["n_open"] for _, p in ok], marker="o", color="steelblue",
                         linewidth=1.4, markersize=8, markeredgecolor="k",
                         markeredgewidth=0.5, label="compromiso re-optimizado")
            axes[1].plot(xs, [p["Z4"] for _, p in ok], marker="o", color="steelblue",
                         linewidth=1.4, markersize=8, markeredgecolor="k",
                         markeredgewidth=0.5, label="Z4 (VAN total)")

        for ax, ref in ((axes[0], ref_n), (axes[1], ref_z4)):
            if ref is not None:
                ax.axhline(ref, color="grey", linestyle=":", linewidth=1.2,
                           label=f"referencia (punto {REFERENCE_POINT})")
            for lv in bad:                              # niveles sin solución factible
                ax.axvline(lv, color="firebrick", linestyle="--", linewidth=1.2, alpha=0.7)
                ax.annotate("ε infactible", (lv, ax.get_ylim()[0]), rotation=90,
                            fontsize=7, color="firebrick", va="bottom", ha="right")

        axes[0].set_xlabel(param)
        axes[0].set_ylabel("nº de plantas abiertas")
        axes[0].set_title(f"{label}: cardinalidad del compromiso")
        axes[0].grid(alpha=0.3)
        axes[0].yaxis.get_major_locator().set_params(integer=True)
        axes[0].legend(fontsize=7)

        axes[1].set_xlabel(param)
        axes[1].set_ylabel("Z4 — VAN total (EUR)")
        axes[1].set_title(f"{label}: VAN del compromiso")
        axes[1].grid(alpha=0.3)
        axes[1].legend(fontsize=7)

        fig.tight_layout()
        path = os.path.join(outdir, f"comparacion_{sid}.png")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        print(f"  figura comparativa → {path}")


# ===========================================================================
# PRE-VUELO
# ===========================================================================
def report_sizes(runs, zone):
    """Tamaño del modelo por geometría distinta (carga de datos, sin solver). El
    escenario de catchment reescala el grafo entero; los otros tres no tocan
    I/J/K, así que basta una carga por geometría."""
    print("\n  Tamaño del modelo por geometría (carga de datos, sin solver):")
    geo_keys, seen = {}, {}
    for rid, cfg in runs.items():
        key = (cfg["D_BAR_IJ"], cfg["D_BAR_JK"], cfg["D_BAR_KK_KM"], cfg["D_BAR_CENTER_K"])
        geo_keys.setdefault(key, []).append(rid)

    for key, rids in sorted(geo_keys.items()):
        cfg = dict(runs[rids[0]])
        try:
            apply_config(cfg, "sizecheck", TIME_LIMIT_DEFAULT, MIP_GAP_DEFAULT)
            D = sens.load_instance(zone)
            if D is None:
                raise RuntimeError("load_instance devolvió None")
            n_b = len(D["A_IJ"]) * len(D["M"]) * len(D["TT"])
            n_x = len(D["A_JK"]) * len(D["M"]) * len(D["TT"])
            gen = sum(D["Gamma"].values())
            seen[key] = n_b
            warn = "   <-- MUY GRANDE" if n_b > 1_500_000 else ""
            print(f"    d̄^IJ={key[0]:<5}  |I|={len(D['I']):>6}  |J|={len(D['J']):>4}  "
                  f"|K|={len(D['K']):>3}  b≈{n_b:>10,}  x≈{n_x:>8,}  "
                  f"gen={gen/1e6:>6.2f}M kg{warn}")
            print(f"          corridas: {', '.join(rids)}")
        except Exception as exc:
            print(f"    d̄^IJ={key[0]:<5}  [ERROR] {type(exc).__name__}: {exc}")
            print(f"          corridas: {', '.join(rids)}")

    if any(v > 1_500_000 for v in seen.values()):
        print("\n  ⚠ Alguna geometría supera ~1.5M variables de flujo b (es el nivel alto de "
              "\n    catchment). Cambiar de zona NO es opción aquí: los ε del punto son "
              "\n    absolutos y sólo tienen sentido en la zona de la referencia. Si esa celda "
              "\n    no cierra, súbele el presupuesto (--time-limit) o córrela sola con "
              "\n    --only catchment y déjala toda la noche; el resto del barrido no depende "
              "\n    de ella y las tablas se escriben tras cada corrida.")


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    global REFERENCE_DIR, REFERENCE_POINT
    ap = argparse.ArgumentParser(
        description="Sensibilidad OFAT de la solución de compromiso (instancia Nancy).")
    ap.add_argument("--zone", default=ZONE_DEFAULT, help=f"Zona a resolver (def. {ZONE_DEFAULT})")
    ap.add_argument("--outdir", default=None, help="Carpeta de salida del barrido")
    ap.add_argument("--reference", default=REFERENCE_DIR,
                    help=f"Corrida con el frente de referencia (def. {REFERENCE_DIR})")
    ap.add_argument("--point", type=int, default=REFERENCE_POINT,
                    help=f"Índice 1-based del punto del frente (def. {REFERENCE_POINT})")
    ap.add_argument("--time-limit", dest="time_limit", type=int, default=TIME_LIMIT_DEFAULT,
                    help=f"Límite de tiempo por solve, s (def. {TIME_LIMIT_DEFAULT})")
    ap.add_argument("--mip-gap", dest="mip_gap", type=float, default=MIP_GAP_DEFAULT,
                    help=f"Gap relativo del solver (def. {MIP_GAP_DEFAULT})")
    ap.add_argument("--augw", type=float, default=sens.ECON_AUGW,
                    help=f"eps del término aumentado, ec. (6) (def. {sens.ECON_AUGW})")
    ap.add_argument("--only", nargs="+", default=None, metavar="ESC",
                    help="Resolver sólo estos escenarios: " + " ".join(s for s, _, _, _ in SWEEPS))
    ap.add_argument("--force", action="store_true",
                    help="Re-resolver corridas que ya tienen compromiso_results.json")
    ap.add_argument("--allow-zone-mismatch", dest="allow_zone_mismatch", action="store_true",
                    help="Permitir una zona distinta a la de la corrida de referencia "
                         "(por defecto se aborta: los ε son absolutos y no trasladan)")
    ap.add_argument("--dry-run", dest="dry_run", action="store_true",
                    help="Imprime el plan y los tamaños de modelo, sin resolver")
    args = ap.parse_args()

    REFERENCE_DIR = args.reference
    REFERENCE_POINT = args.point

    runs, plan = build_plan()
    if args.only:
        bad = [s for s in args.only if s not in {x[0] for x in SWEEPS}]
        if bad:
            raise SystemExit(f"Escenario(s) desconocido(s): {bad}. "
                             f"Opciones: {[s for s, _, _, _ in SWEEPS]}")
        plan = [r for r in plan if r[0] in set(args.only) | {"base"}]
        keep = {r[4] for r in plan}
        runs = {rid: cfg for rid, cfg in runs.items() if rid in keep}

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    outdir = args.outdir or os.path.join(
        sens.REPO_ROOT, "results",
        f"results_compromiso_{args.zone}_{ts}")

    print(f"\n=== Sensibilidad OFAT de la solución de compromiso | zona: {args.zone} ===")
    print(f"  outdir     : {outdir}")
    print(f"  método     : se re-resuelve el subproblema AUGMECON del punto, con los ε "
          f"FIJOS en valor absoluto")
    print(f"  solver     : TL={args.time_limit}s  gap={args.mip_gap}  augw={args.augw}")
    print(f"  corridas   : {len(runs)} configuraciones × 1 solve = {len(runs)} solves")
    print()

    ref_pt, eps, ranges, ref_params = load_reference_point(REFERENCE_DIR, REFERENCE_POINT)

    # Los ε son valores ABSOLUTOS calibrados sobre el frente de una zona concreta.
    # En una zona más chica, Z3 sencillamente no llega a esa escala (p. ej. el
    # máximo alcanzable en `centro` es ~17× menor que el eps_Z3 del punto 21), así
    # que TODAS las configuraciones darían infactible — no por el parámetro
    # perturbado, sino por un desajuste de escala. Se aborta antes de quemar horas.
    ref_zone = ref_params.get("ZONE")
    if ref_zone and ref_zone != args.zone and not args.allow_zone_mismatch:
        raise SystemExit(
            f"\nLa corrida de referencia es de la zona '{ref_zone}' y estás pidiendo "
            f"'{args.zone}'.\nLos ε del punto {REFERENCE_POINT} "
            f"({'  '.join(f'{k}={v:,.1f}' for k, v in eps.items())}) son valores absolutos "
            f"de la escala de '{ref_zone}': en otra zona el frente vive en otra escala y "
            f"las 9 configuraciones saldrían infactibles por construcción, no por el "
            f"parámetro barrido.\n\nOpciones:\n"
            f"  · correr en la zona de la referencia:  --zone {ref_zone}\n"
            f"  · usar un frente de referencia de '{args.zone}':  --reference <carpeta>\n"
            f"  · forzarlo igual (sabiendo lo anterior):  --allow-zone-mismatch\n")
    if ref_zone and ref_zone != args.zone:
        print(f"\n  ⚠ zona '{args.zone}' ≠ zona de la referencia '{ref_zone}' "
              f"(--allow-zone-mismatch): es muy probable que los ε sean inalcanzables.")

    print()
    for sid, label, param, lv, run_id, is_base in plan:
        mark = "  (= base, se reutiliza)" if is_base and sid != "base" else ""
        lvtxt = "—" if lv is None else f"{lv}"
        print(f"    {sid:10s} {param:24s} {lvtxt:>6s}  ->  {run_id}{mark}")

    check_instance_cache(args.zone, list(runs.values()))

    if args.dry_run:
        report_sizes(runs, args.zone)
        print("\n--dry-run: no se resuelve nada. Quita la bandera para lanzar el barrido.\n")
        return

    runs_dir = os.path.join(outdir, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    results = {}
    t_all = time.time()

    for n, (run_id, cfg) in enumerate(runs.items(), 1):
        rundir = os.path.join(runs_dir, run_id)
        done = os.path.exists(os.path.join(rundir, "compromiso_results.json"))

        print(f"\n{'='*78}\n[{n}/{len(runs)}] configuración '{run_id}'"
              f"{'  (ya resuelta: se reutiliza)' if done and not args.force else ''}\n{'='*78}")
        diff = "  ".join(f"{k}={v}" for k, v in cfg.items() if BASE_CFG.get(k) != v)
        print(f"  config: {diff if diff else 'caso base (sin desviaciones)'}")

        if done and not args.force:
            status = dict(run_id=run_id, status="reutilizada", error=None, wall_time_s="")
        else:
            status = solve_run(run_id, cfg, args.zone, rundir, ref_pt, eps, ranges,
                               args.time_limit, args.mip_gap, args.augw, REFERENCE_POINT)
            print(f"  -> {status['status']}  ({status.get('wall_time_s', '?')}s)")

        loaded = load_run(rundir) or {}
        results[run_id] = dict(status=status,
                               solucion=loaded.get("solucion"),
                               factible=loaded.get("factible", False),
                               diagnostico_infactibilidad=loaded.get("diagnostico_infactibilidad"),
                               geo=loaded.get("geo", {}))

        write_aggregates(outdir, plan, results, runs, ref_pt, eps, REFERENCE_POINT)

    make_comparison_plots(outdir, plan, results, ref_pt)
    write_aggregates(outdir, plan, results, runs, ref_pt, eps, REFERENCE_POINT)

    fact = [rid for rid, r in results.items() if r.get("factible")]
    infact = [rid for rid, r in results.items()
              if r["status"].get("status") == "infactible"]
    err = [rid for rid, r in results.items() if r["status"].get("status") == "error"]
    print(f"\n{'='*78}")
    print(f"Barrido terminado en {(time.time()-t_all)/60:.1f} min — "
          f"{len(fact)}/{len(runs)} configuraciones con solución factible.")
    if infact:
        print(f"  ε infactible (la especificación del compromiso no se sostiene): {infact}")
        print(f"    -> ver diagnostico_infactibles.csv para saber QUÉ ε no se alcanza")
    if err:
        print(f"  con error: {err}  (ver runs/<id>/run.log)")
    print(f"\nTablas agregadas en {outdir}:")
    for f in ("sensibilidad_compromiso.csv", "cardinalidad_compromiso.csv",
              "diagnostico_infactibles.csv", "plantas_detalle.csv", "runs_manifest.json"):
        print(f"    {f}")
    print()


if __name__ == "__main__":
    main()
