# -*- coding: utf-8 -*-
"""
sensibilidad_articulo.py — análisis de sensibilidad de un factor (OFAT) sobre las
SOLUCIONES MONO-OBJETIVO de la instancia Nancy.

Qué hace y qué NO hace
----------------------
NO reconstruye el frente de Pareto en cada nivel de parámetro (eso serían ~450
solves y no es lo que pide un análisis de sensibilidad). Toma las tres soluciones
mono-objetivo de la corrida de referencia y, para cada nivel de cada parámetro,
LAS VUELVE A OPTIMIZAR — que es la definición estándar de sensibilidad de un
factor (ceteris paribus): se perturba un parámetro, se re-resuelve, y se observa
cómo se mueve el óptimo. Re-optimizar (y no fijar la solución y re-evaluarla) es
lo que permite que cambie la cardinalidad, que es justo lo que el manuscrito
quiere someter a robustez ("the number and identity of facilities that open").

Comparabilidad con la referencia
--------------------------------
Los anclas de la corrida de referencia NO son óptimos mono-objetivo ingenuos:
son óptimos LEXICOGRÁFICOS (Mavrotas 2009, Sec. 3.1), construidos para no quedar
en un óptimo alternativo arbitrario del MILP. Este script re-optimiza con LA
MISMA cadena lexicográfica, así que cada solución perturbada es directamente
comparable con su ancla de referencia:

    fila Z4 (max VAN total) :  Z4 -> Z2 -> Z3
    fila Z2 (min GEI)       :  Z2 -> Z4 -> Z3
    fila Z3 (max cobertura) :  Z3 -> Z4 -> Z2

Con --no-lex se hace un solve mono-objetivo puro (1 solve por objetivo en vez de
3). Es más barato pero las soluciones dejan de ser comparables con la referencia
y quedan expuestas a óptimos alternativos degenerados; usarlo sólo para tantear.

Referencia
----------
    results/pareto_front/econ_results.json   (zona full, escenario BASE)

      ancla            objetivo          n_open  plantas                   gap
      lex_anchor_Z4    max VAN total          3  [2, 7, 8]              42.80%  <-- ojo
      lex_anchor_Z2    min GEI                8  [1,4,7,8,9,10,12,14]    1.84%
      lex_anchor_Z3    max cobertura          8  [1,4,7,8,9,11,12,14]    0.00%

⚠ El ancla Z4 de la referencia cerró con 42.8% de gap (topó el límite de 1000 s).
El "3 plantas" del óptimo económico NO está probado: es un incumbente, no un
óptimo. Cualquier lectura de robustez de la cardinalidad sobre esa fila hereda
esa incertidumbre. Las columnas gap_pct de las salidas lo reportan en cada nivel.

Cada solución de referencia se usa además como WARM START de su fila, que es la
otra razón por la que la sensibilidad se ancla "sobre las soluciones mono-obj".

Escenarios (tabla del manuscrito)
---------------------------------
  Escenario                   Parámetro                     Base    Niveles
  --------------------------  ----------------------------  ------  ---------------------
  Price volatility            PRICE_FACTOR   (p_m, todos)   1.0     ×{0.7, 1.0, 1.3}
  Catchment sensitivity       D_BAR_IJ       (d̄^IJ)         0.25    {0.10, 0.25, 0.50} km
  Municipal credit intensity  MUNICIPAL_CREDIT_SHARE         0.25    {0.00, 0.25, 0.50}
  Facility cost/capacity      F_fix y Q_prod (F_k, Q̄^prod)  1.0     ×{0.7, 1.0, 1.3}

El nivel que coincide con el base no se re-resuelve: las cuatro filas base
apuntan a la única corrida `base`. Son 9 configuraciones, no 13, y por cada una
se resuelven 3 filas lexicográficas => 9 × 3 × 3 ≈ 81 solves (contra ~450 del
barrido multi-objetivo completo).

USO
---
  python sensibilidad_articulo.py --dry-run       # plan + tamaño de modelo, sin resolver
  python sensibilidad_articulo.py                 # barrido completo, zona full
  python sensibilidad_articulo.py --zone centro   # prueba rápida en zona chica
  python sensibilidad_articulo.py --only price catchment
  python sensibilidad_articulo.py --objectives Z4          # sólo la fila económica
  python sensibilidad_articulo.py --time-limit 900 --force

Reanuda por defecto: una configuración cuyo `monoobj_results.json` ya existe se
salta y se relee. `--force` la re-resuelve.

SALIDAS (bajo --outdir, por defecto results_articulo_<zona>_<timestamp>/)
--------------------------------------------------------------------------
  runs/<run_id>/monoobj_results.json   las 3 soluciones completas + config + payoff
  runs/<run_id>/detalle_csv/*.csv      plantas, hubs, flujos, zonas (por solución)
  runs/<run_id>/run_config.json        configuración exacta de esa corrida
  runs/<run_id>/run.log                stdout completo del solve
  sensibilidad_monoobj.csv    1 fila por (escenario, nivel, objetivo) — TABLA PRINCIPAL
  cardinalidad.csv            foco en n_open y en qué plantas entran/salen vs referencia
  plantas_frecuencia.csv      1 fila por (escenario, nivel, objetivo, planta)
  payoff_tables.csv           tabla de pagos de cada configuración
  runs_manifest.json          config + estado + tiempos de cada corrida
  comparacion_*.png           cardinalidad y valor objetivo por nivel

Las tablas agregadas se reescriben DESPUÉS DE CADA CORRIDA: un barrido
interrumpido deja en disco tablas válidas con lo ya resuelto.
"""

import os
import sys
import csv
import json
import time
import argparse
import datetime
import traceback
from collections import Counter

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pyomo.environ as pyo
import sensibilidad as sens


# ===========================================================================
# CONFIGURACIÓN DEL BARRIDO
# ===========================================================================

# Corrida de referencia: de aquí salen las soluciones mono-objetivo que se
# sensibilizan (y que sirven de warm start y de término de comparación).
REFERENCE_DIR = os.path.join(sens.REPO_ROOT, "results", "pareto_front")

# Objetivos mono-objetivo a sensibilizar = las tres filas de la tabla de pagos
# de la referencia. Z1 (VAN max-min) y Z5 (output) quedan como diagnóstico, igual
# que en la corrida de referencia.
OBJECTIVES = ["Z4", "Z2", "Z3"]
PAYOFF_OBJS = ("Z4", "Z2", "Z3")     # objetivos que participan de la cadena lexicográfica

# Límite de tiempo por solve. Al no haber grilla ε son ~81 solves en vez de ~450,
# así que se puede ser más generoso que en un barrido multi-objetivo. El warm
# start desde la solución de referencia ayuda mucho en los niveles cercanos al base.
TIME_LIMIT_DEFAULT = 600
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
    FACILITY_FACTOR=1.0,      # escala F_k y Q̄^prod_k CONJUNTAMENTE (ver scale_facilities)
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


def _g(p, key, default=""):
    """Lectura tolerante de un punto que pudo pasar por JSON (claves str/int)."""
    v = p.get(key, default)
    return v


def check_instance_cache(zone, cfg_list):
    """El caché de instancias sólo guarda arcos hasta arc_ij_max / arc_jk_max
    (build_instances.py). Un D_BAR_IJ por encima de ese techo NO amplía el grafo:
    devolvería el mismo grafo recortado y el barrido de catchment saldría plano
    por un artefacto del caché, no por el modelo."""
    man_path = os.path.join(sens.INSTANCES_DIR, "manifest.json")
    if not os.path.exists(man_path):
        print(f"  [warning] Sin {man_path}: no puedo verificar el techo de arcos del caché.")
        return
    with open(man_path, encoding="utf-8") as fh:
        zones = json.load(fh).get("zones", [])
    z = next((x for x in zones if x.get("name") == zone), None)
    if z is None:
        print(f"  [warning] La zona '{zone}' no está en el manifest; no verifico el techo de arcos.")
        return

    ij_max = float(z.get("arc_ij_max", 0.0))
    jk_max = float(z.get("arc_jk_max", 0.0))
    need_ij = max(c["D_BAR_IJ"] for c in cfg_list)
    need_jk = max(c["D_BAR_JK"] for c in cfg_list)

    if need_ij > ij_max + 1e-9:
        raise SystemExit(
            f"\nD_BAR_IJ máximo del barrido = {need_ij} km, pero el caché de '{zone}' sólo "
            f"guarda arcos I-J hasta {ij_max} km.\nRegenera el caché (build_instances.py, "
            f"ajustando arc_ij_max) — si no, el nivel alto de 'catchment' sería idéntico "
            f"al techo del caché y el escenario saldría plano por artefacto.")
    if need_jk > jk_max + 1e-9:
        raise SystemExit(
            f"\nD_BAR_JK = {need_jk} km supera el techo del caché ({jk_max} km). "
            f"Regenera con build_instances.py.")

    if abs(need_ij - ij_max) < 1e-9:
        print(f"  [info] El nivel alto de catchment ({need_ij} km) coincide EXACTAMENTE con el "
              f"techo del caché (arc_ij_max={ij_max}). Es válido, pero no se puede subir más "
              f"sin re-correr build_instances.py.")


# ===========================================================================
# REFERENCIA
# ===========================================================================
def load_reference(refdir):
    """Lee las soluciones mono-objetivo (anclas lexicográficas) de la corrida de
    referencia. Devuelve (anchors, params) con anchors = {objetivo: punto}."""
    path = os.path.join(refdir, "econ_results.json")
    if not os.path.isabs(path):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
    if not os.path.exists(path):
        raise SystemExit(f"\nNo existe la corrida de referencia: {path}\n"
                         f"Ajusta REFERENCE_DIR o pasa --reference <carpeta>.")
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)

    anchors = {}
    for p in data.get("raw", []):
        st = str(p.get("status", ""))
        if st.startswith("lex_anchor_"):
            anchors[st.replace("lex_anchor_", "")] = p

    if not anchors:
        raise SystemExit(f"\n{path} no contiene anclas lex_anchor_*. ¿Es una corrida "
                         f"anterior al payoff lexicográfico?")

    print(f"  referencia: {path}")
    for obj in OBJECTIVES:
        a = anchors.get(obj)
        if a is None:
            print(f"    {obj}: (ausente en la referencia)")
            continue
        gp = a.get("gap_pct")
        gtxt = f"{gp:.2f}%" if isinstance(gp, (int, float)) else str(gp)
        warn = "   <-- gap alto: el óptimo NO está probado" if isinstance(gp, (int, float)) and gp > 5 else ""
        print(f"    {obj}: n_open={a.get('n_open')}  plantas={a.get('opened')}  "
              f"{obj}={a.get(obj, 0):,.1f}  gap={gtxt}{warn}")
    return anchors, data.get("params", {})


def seed_from_reference(m, D, ref_pt):
    """Siembra y[k,t] con la solución de referencia, como MIP start.

    Es la conexión concreta con "la sensibilidad se hace sobre las soluciones
    mono-obj": cada fila arranca desde su propia solución de referencia. Sólo se
    siembran plantas que EXISTEN en la instancia perturbada — al mover d̄^IJ el
    conjunto K cambia (9/10/10 plantas según el nivel), así que una planta de la
    referencia puede no estar presente.

    Es un arranque PARCIAL (sólo las binarias): Gurobi completa el resto. Si el
    backend no lo acepta, `_safe_solve` reintenta solo, sin warmstart."""
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
        if k not in K:                                 # planta ausente de esta geometría
            continue
        for t_str, v in serie.items():
            t = int(t_str)
            if t not in TT:                            # horizonte distinto al de la referencia
                continue
            m.y[k, t].value = int(v)
            seeded += 1
    return seeded


# ===========================================================================
# APLICACIÓN DE LA CONFIGURACIÓN
# ===========================================================================
def apply_config(cfg, run_id, time_limit, mip_gap):
    """Escribe la configuración de la corrida en los globals de sensibilidad.py.

    Funciona porque `load_instance`, `_build_final_model_dict` y `build_model`
    leen esos globals EN TIEMPO DE LLAMADA, no al importar. Lo único sin global
    es el par (F_k, Q̄^prod), que se escala sobre D — ver scale_facilities."""
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


def scale_facilities(D, factor):
    """Escenario 4: escala F_k y Q̄^prod_k CONJUNTAMENTE, como pide el manuscrito.

    No hay global para esto: F_fix y Q_prod se leen directo del caché en
    `_build_final_model_dict`. Se escalan sobre el dict D ya construido, que es
    lo que consume `build_model`:
        · D["Q_prod"][k] -> restricción intake_cap (Σx ≤ Q_prod·y)
        · D["F_fix"][k]  -> costo fijo en π_k
    Son sus dos únicos consumidores, así que equivale a escalar la fuente.

    M_pi_global (big-M de Z1) se recalcula porque depende de max(Q_prod): si se
    escala Q_prod hacia arriba sin tocarlo, el big-M queda POR DEBAJO del π_k
    alcanzable y phi_constraint recortaría soluciones factibles."""
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
# RE-OPTIMIZACIÓN MONO-OBJETIVO
# ===========================================================================
def solve_mono(m, solver, primary, ref_pt, D, use_lex=True):
    """Re-optimiza UNA solución mono-objetivo bajo la configuración vigente.

    Con use_lex=True replica la fila `primary` de la tabla de pagos lexicográfica
    de Mavrotas (2009), Sec. 3.1 — la MISMA construcción que produjo el ancla de
    referencia, de modo que ambas son comparables:
        1. optimiza `primary`;
        2. amarra ese óptimo con una restricción de holgura atada al gap del
           solver (`_add_lex_fix`) y optimiza el siguiente objetivo del orden;
        3. repite; el punto final de la cadena es el ancla.
    Con use_lex=False hace un solo solve de `primary` (más barato, pero expuesto
    a óptimos alternativos: dos niveles vecinos pueden devolver redes distintas
    con idéntico valor objetivo, y la comparación de cardinalidad se vuelve ruido).

    Devuelve (punto, meta) o (None, meta) si no hubo solución."""
    order = sens._lex_row_order(primary, PAYOFF_OBJS) if use_lex else [primary]
    n_seed = seed_from_reference(m, D, ref_pt)
    fixes = []
    z_prev = None
    meta = {"gap_pct": "NA", "bound": None, "wall_time_s": 0.0}
    cat = "nosolution"

    print(f"    cadena {'→'.join(order)}"
          + (f"  (warm start: {n_seed} binarias desde la referencia)" if n_seed else ""))

    for step, key in enumerate(order):
        if step > 0:
            name = f"lex_fix_{primary}_{step}"
            sens._add_lex_fix(m, name, order[step - 1], z_prev)
            fixes.append(name)
        sens.set_objective(m, sens._obj_expr(m, key), sens.OBJ_SENSE[key])
        # warmstart en el paso 0 = la solución de referencia sembrada arriba;
        # en los pasos siguientes = el incumbente del paso anterior de la cadena.
        cat, _, meta = sens._safe_solve(solver, m, tee=sens.SOLVER_TEE, warmstart=True)
        gp = meta["gap_pct"]
        gtxt = f"{gp:.2f}%" if isinstance(gp, (int, float)) else str(gp)
        if cat in ("infeasible", "nosolution", "error"):
            print(f"      paso {step} ({key}): {cat} — se corta la cadena y se "
                  f"registra el último punto factible")
            break
        z_prev = float(pyo.value(sens._obj_expr(m, key)))
        print(f"      paso {step} ({key}): {cat:<9} {key}={z_prev:>14,.1f}  "
              f"gap={gtxt:<8} t={meta['wall_time_s']:>5.0f}s")

    pt = None
    try:
        pt = sens._record_point(m, status=f"mono_{primary}",
                                gap_pct=meta["gap_pct"], bound=meta["bound"],
                                wall_time_s=meta["wall_time_s"],
                                objective=primary, lex_order=order,
                                termination=cat)
    except Exception as exc:
        print(f"      [warning] no se pudo registrar el punto: {exc}")

    for name in fixes:                       # la fila siguiente parte limpia
        m.del_component(name)
    return pt, meta


# ===========================================================================
# PLAN DE CORRIDAS
# ===========================================================================
def build_plan():
    """Devuelve (runs, plan).

    runs: run_id -> cfg   — configuraciones ÚNICAS a resolver.
    plan: (scenario_id, etiqueta, parámetro, nivel, run_id, es_base) — las filas
          de la tabla del manuscrito; varias apuntan al mismo run_id porque el
          nivel base se comparte entre los cuatro escenarios."""
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
# UNA CORRIDA (= una configuración de parámetros, 3 filas mono-objetivo)
# ===========================================================================
def solve_run(run_id, cfg, zone, rundir, ref_anchors, objectives,
              time_limit, mip_gap, use_lex):
    """Resuelve las filas mono-objetivo bajo una configuración. Nunca lanza: los
    errores se capturan para que una corrida caída no aborte el barrido."""
    os.makedirs(rundir, exist_ok=True)

    cfg_payload = dict(
        run_id=run_id, zone=zone, config=cfg,
        reference=REFERENCE_DIR, objectives=objectives, lexicographic=use_lex,
        lex_chains={o: sens._lex_row_order(o, PAYOFF_OBJS) for o in objectives} if use_lex else None,
        solver=dict(time_limit_s=time_limit, mip_gap=mip_gap),
        horizon=sens.HORIZON, alpha=sens.ALPHA,
        yield_mode=sens.YIELD_MODE, use_disposal=sens.USE_DISPOSAL,
        lca_credit_m=sens.LCA_CREDIT_M, c_mun=sens.C_MUN,
        started_at=datetime.datetime.now().isoformat(timespec="seconds"),
    )
    with open(os.path.join(rundir, "run_config.json"), "w", encoding="utf-8") as fh:
        json.dump(cfg_payload, fh, indent=1, ensure_ascii=False)

    tee = _Tee(sys.stdout, os.path.join(rundir, "run.log"))
    real_stdout = sys.stdout
    sys.stdout = tee
    t0 = time.time()
    status = dict(run_id=run_id, status="ok", error=None)
    points = {}

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

        for obj in objectives:
            print(f"\n  == fila mono-objetivo {obj} "
                  f"({sens.OBJ_LABEL.get(obj, '')}, {'max' if sens.OBJ_SENSE[obj] == pyo.maximize else 'min'}) ==")
            pt, _ = solve_mono(m, solver, obj, ref_anchors.get(obj), D, use_lex=use_lex)
            if pt is None:
                print(f"    -> sin solución para {obj}")
                continue
            points[obj] = pt
            print(f"    -> n_open={pt['n_open']}  plantas={pt['opened']}  "
                  f"Z4={pt.get('Z4', 0):,.0f}  Z2={pt.get('Z2', 0):,.0f}  Z3={pt.get('Z3', 0):,.1f}")

        # Tabla de pagos de ESTA configuración: una fila por objetivo, con el
        # valor de todos los objetivos en el punto final de su cadena.
        payoff = {o: {k: float(p.get(k, 0.0)) for k in PAYOFF_OBJS}
                  for o, p in points.items()}

        save_run(points, payoff, D, cfg, run_id, zone, rundir, objectives, use_lex)
        status.update(n_soluciones=len(points))

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


def save_run(points, payoff, D, cfg, run_id, zone, rundir, objectives, use_lex):
    """Guarda las soluciones de una configuración: JSON completo + los CSV
    detallados de sensibilidad.py (que operan sobre puntos de `_record_point`,
    así que sirven igual para soluciones mono-objetivo)."""
    ordered = [points[o] for o in objectives if o in points]
    payload = dict(
        run_id=run_id, zone=zone, config=cfg, reference=REFERENCE_DIR,
        lexicographic=use_lex, objectives=objectives,
        payoff=payoff,
        soluciones={o: sens._strip_heavy(p, keep_full=True) for o, p in points.items()},
        geo=sens._build_geo(D),
    )
    path = os.path.join(rundir, "monoobj_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False)
    print(f"\n  resultados → {path}")

    if ordered:
        try:
            sens.save_detailed_csv(ordered, D, rundir)
        except Exception as exc:
            print(f"  [warning] save_detailed_csv falló: {exc}")


def load_run(rundir):
    """Relee una corrida ya resuelta (modo reanudar)."""
    path = os.path.join(rundir, "monoobj_results.json")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


# ===========================================================================
# TABLAS AGREGADAS
# ===========================================================================
def _set_delta(ref_opened, new_opened):
    """Qué plantas entran y salen respecto de la solución de referencia."""
    r, n = set(int(x) for x in ref_opened or []), set(int(x) for x in new_opened or [])
    entran = sorted(n - r)
    salen = sorted(r - n)
    jac = (len(r & n) / len(r | n)) if (r | n) else 1.0
    return entran, salen, round(jac, 4)


def write_aggregates(outdir, plan, results, runs, ref_anchors, objectives):
    """Reescribe las tablas agregadas + el manifest. Se llama tras CADA corrida."""

    def _w(name, header, rows):
        with open(os.path.join(outdir, name), "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(header)
            wr.writerows(rows)

    # --- 1. TABLA PRINCIPAL: una fila por (escenario, nivel, objetivo) --------
    hdr = ["escenario", "etiqueta", "parametro", "nivel", "es_base", "run_id", "estado",
           "objetivo", "sentido", "valor_objetivo",
           "Z1_VAN_peor_EUR", "Z2_GHG_kgCO2e", "Z3_cobertura", "Z4_NPV_EUR", "Z5_output_kg",
           "cobertura_pct", "colectado_total_kg",
           "n_plantas", "plantas",
           "ref_n_plantas", "ref_plantas", "plantas_entran", "plantas_salen", "jaccard_vs_ref",
           "delta_objetivo_vs_ref", "delta_pct_vs_ref",
           "gap_pct", "bound", "wall_time_s", "terminacion", "error"]
    rows = []
    for sid, label, param, lv, run_id, is_base in plan:
        r = results.get(run_id, {})
        st = r.get("status", {})
        sols = r.get("soluciones", {})
        for obj in objectives:
            p = sols.get(obj)
            ref = ref_anchors.get(obj, {})
            if p is None:
                rows.append([sid, label, param, ("" if lv is None else lv), int(is_base),
                             run_id, st.get("status", "pendiente"), obj,
                             "max" if sens.OBJ_SENSE[obj] == pyo.maximize else "min"]
                            + [""] * 19 + [st.get("error", "") or ""])
                continue
            val = float(p.get(obj, 0.0))
            rval = float(ref.get(obj, 0.0)) if ref else None
            entran, salen, jac = _set_delta(ref.get("opened"), p.get("opened"))
            dlt = (val - rval) if rval is not None else ""
            dpct = (100.0 * (val - rval) / abs(rval)) if rval not in (None, 0) else ""
            rows.append([
                sid, label, param, ("" if lv is None else lv), int(is_base), run_id,
                st.get("status", "pendiente"), obj,
                "max" if sens.OBJ_SENSE[obj] == pyo.maximize else "min", val,
                p.get("Z1", ""), p.get("Z2", ""), p.get("Z3", ""),
                p.get("Z4", p.get("NPV", "")), p.get("Z5", ""),
                p.get("coverage_pct", ""), p.get("collected_total", ""),
                p.get("n_open", ""), " ".join(map(str, p.get("opened", []))),
                ref.get("n_open", ""), " ".join(map(str, ref.get("opened", []))),
                " ".join(map(str, entran)), " ".join(map(str, salen)), jac,
                dlt, dpct,
                p.get("gap_pct", "NA"), p.get("bound", ""),
                p.get("wall_time_s", ""), p.get("termination", ""), "",
            ])
    _w("sensibilidad_monoobj.csv", hdr, rows)

    # --- 2. CARDINALIDAD: el resultado estructural del manuscrito -------------
    hdr = ["escenario", "parametro", "nivel", "objetivo",
           "n_plantas", "ref_n_plantas", "delta_n_plantas",
           "plantas", "plantas_entran", "plantas_salen", "jaccard_vs_ref",
           "misma_cardinalidad", "misma_identidad", "gap_pct"]
    rows = []
    for sid, label, param, lv, run_id, is_base in plan:
        sols = results.get(run_id, {}).get("soluciones", {})
        for obj in objectives:
            p = sols.get(obj)
            if p is None:
                continue
            ref = ref_anchors.get(obj, {})
            entran, salen, jac = _set_delta(ref.get("opened"), p.get("opened"))
            n, rn = p.get("n_open"), ref.get("n_open")
            rows.append([sid, param, ("" if lv is None else lv), obj,
                         n, rn, (n - rn) if isinstance(n, int) and isinstance(rn, int) else "",
                         " ".join(map(str, p.get("opened", []))),
                         " ".join(map(str, entran)), " ".join(map(str, salen)), jac,
                         int(n == rn) if isinstance(n, int) and isinstance(rn, int) else "",
                         int(not entran and not salen),
                         p.get("gap_pct", "NA")])
    _w("cardinalidad.csv", hdr, rows)

    # --- 3. frecuencia de apertura por planta ---------------------------------
    hdr = ["escenario", "nivel", "run_id", "objetivo", "k", "nombre", "x", "y",
           "abierta", "apertura_t", "VAN_planta_EUR", "en_referencia"]
    rows = []
    for sid, label, param, lv, run_id, is_base in plan:
        r = results.get(run_id, {})
        coords = r.get("geo", {}).get("coords_k", {})
        for obj in objectives:
            p = r.get("soluciones", {}).get(obj)
            if p is None:
                continue
            refset = set(int(x) for x in (ref_anchors.get(obj, {}).get("opened") or []))
            pi = p.get("pi_individual", {}) or {}
            for k in p.get("opened", []):
                ck = coords.get(str(k), [None, None, ""])
                rows.append([sid, ("" if lv is None else lv), run_id, obj, k,
                             (ck[2] if len(ck) > 2 else ""), ck[0], ck[1], 1,
                             p.get("open_period", {}).get(str(k),
                                   p.get("open_period", {}).get(k, "")),
                             pi.get(str(k), pi.get(k, "")),
                             int(int(k) in refset)])
    _w("plantas_frecuencia.csv", hdr, rows)

    # --- 4. tablas de pagos por configuración ---------------------------------
    hdr = ["escenario", "nivel", "run_id", "fila_objetivo", "objetivo_evaluado", "valor"]
    rows = []
    for sid, label, param, lv, run_id, is_base in plan:
        for prim, cellrow in (results.get(run_id, {}).get("payoff", {}) or {}).items():
            if not isinstance(cellrow, dict):
                continue
            for o, v in cellrow.items():
                rows.append([sid, ("" if lv is None else lv), run_id, prim, o, v])
    _w("payoff_tables.csv", hdr, rows)

    # --- 5. manifest -----------------------------------------------------------
    manifest = dict(
        generated_at=datetime.datetime.now().isoformat(timespec="seconds"),
        analysis="one-factor (OFAT) sensitivity, re-optimizing the mono-objective "
                 "(lexicographic) solutions of the reference run",
        reference=REFERENCE_DIR,
        reference_anchors={o: dict(n_open=a.get("n_open"), opened=a.get("opened"),
                                   value=a.get(o), gap_pct=a.get("gap_pct"))
                           for o, a in ref_anchors.items()},
        base_config=BASE_CFG,
        sweeps=[dict(id=s, label=l, param=p, levels=v) for s, l, p, v in SWEEPS],
        objectives=objectives,
        lex_chains={o: sens._lex_row_order(o, PAYOFF_OBJS) for o in objectives},
        runs={rid: dict(config=cfg, status=results.get(rid, {}).get("status", {}))
              for rid, cfg in runs.items()},
        plan=[dict(escenario=s, etiqueta=l, parametro=p, nivel=lv, run_id=r, es_base=b)
              for s, l, p, lv, r, b in plan],
    )
    with open(os.path.join(outdir, "runs_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1, ensure_ascii=False)


# ===========================================================================
# GRÁFICOS COMPARATIVOS
# ===========================================================================
def make_comparison_plots(outdir, plan, results, ref_anchors, objectives):
    """Por escenario: cardinalidad vs nivel y valor objetivo vs nivel, una serie
    por objetivo. La línea SÍ tiene sentido aquí (es un barrido de un parámetro
    continuo, no un frente de Pareto), pero se marcan los puntos resueltos."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print("  matplotlib no disponible; omito los gráficos comparativos.")
        return

    colors = {"Z4": "steelblue", "Z2": "seagreen", "Z3": "indianred",
              "Z1": "slategrey", "Z5": "darkorange"}

    for sid, label, param, levels in SWEEPS:
        entries = [(lv, results.get(rid, {}).get("soluciones", {}))
                   for s, l, p, lv, rid, b in plan if s == sid]
        if not any(sols for _, sols in entries):
            continue

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        for obj in objectives:
            xs = [lv for lv, sols in entries if sols.get(obj)]
            ns = [sols[obj].get("n_open") for lv, sols in entries if sols.get(obj)]
            vs = [sols[obj].get(obj) for lv, sols in entries if sols.get(obj)]
            if not xs:
                continue
            c = colors.get(obj, "grey")
            axes[0].plot(xs, ns, marker="o", color=c, linewidth=1.3,
                         markersize=7, markeredgecolor="k", markeredgewidth=0.5,
                         label=f"{obj} ({sens.OBJ_LABEL.get(obj,'')})")
            axes[1].plot(xs, vs, marker="o", color=c, linewidth=1.3,
                         markersize=7, markeredgecolor="k", markeredgewidth=0.5,
                         label=obj)
            ra = ref_anchors.get(obj, {})
            if ra.get("n_open") is not None:
                axes[0].axhline(ra["n_open"], color=c, linestyle=":", linewidth=1, alpha=0.55)

        axes[0].set_xlabel(param)
        axes[0].set_ylabel("nº de plantas abiertas")
        axes[0].set_title(f"{label}: cardinalidad por nivel\n(línea punteada = referencia)")
        axes[0].grid(alpha=0.3)
        axes[0].yaxis.get_major_locator().set_params(integer=True)
        axes[0].legend(fontsize=7)

        axes[1].set_xlabel(param)
        axes[1].set_ylabel("valor del objetivo optimizado")
        axes[1].set_title(f"{label}: valor óptimo por nivel")
        axes[1].grid(alpha=0.3)
        axes[1].legend(fontsize=7)

        fig.tight_layout()
        path = os.path.join(outdir, f"comparacion_{sid}.png")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        print(f"  figura comparativa → {path}")


# ===========================================================================
# PRE-VUELO: TAMAÑO DE CADA CORRIDA
# ===========================================================================
def report_sizes(runs, zone):
    """Carga (sólo datos, sin solver) cada geometría distinta y reporta el tamaño
    del modelo. Vale la pena porque el escenario de catchment NO es una
    perturbación suave: d̄^IJ reescala el grafo entero. Los otros tres no tocan
    I/J/K, así que basta una carga por geometría distinta."""
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
            n_bin = len(D["K"]) * len(D["TT"])
            gen = sum(D["Gamma"].values())
            seen[key] = n_b
            warn = "   <-- MUY GRANDE" if n_b > 1_500_000 else ""
            print(f"    d̄^IJ={key[0]:<5}  |I|={len(D['I']):>6}  |J|={len(D['J']):>4}  "
                  f"|K|={len(D['K']):>3}  b≈{n_b:>10,}  x≈{n_x:>8,}  bin={n_bin:>5,}  "
                  f"gen={gen/1e6:>6.2f}M kg{warn}")
            print(f"          corridas: {', '.join(rids)}")
        except Exception as exc:
            print(f"    d̄^IJ={key[0]:<5}  [ERROR] {type(exc).__name__}: {exc}")
            print(f"          corridas: {', '.join(rids)}")

    if any(v > 1_500_000 for v in seen.values()):
        print("\n  ⚠ Alguna geometría supera ~1.5M variables de flujo b. Sobre `full` eso suele "
              "\n    ser intratable con el presupuesto de tiempo por solve. Considera correr el "
              "\n    barrido por zona más chica (--zone centro | sector_N | sector_NE | ...).")


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    global REFERENCE_DIR
    ap = argparse.ArgumentParser(
        description="Sensibilidad OFAT sobre las soluciones mono-objetivo (instancia Nancy).")
    ap.add_argument("--zone", default=ZONE_DEFAULT, help=f"Zona a resolver (def. {ZONE_DEFAULT})")
    ap.add_argument("--outdir", default=None, help="Carpeta de salida del barrido")
    ap.add_argument("--reference", default=REFERENCE_DIR,
                    help=f"Corrida de referencia con los anclas mono-obj (def. {REFERENCE_DIR})")
    ap.add_argument("--objectives", nargs="+", default=OBJECTIVES, metavar="Z",
                    help=f"Objetivos a sensibilizar (def. {' '.join(OBJECTIVES)})")
    ap.add_argument("--time-limit", dest="time_limit", type=int, default=TIME_LIMIT_DEFAULT,
                    help=f"Límite de tiempo por solve, s (def. {TIME_LIMIT_DEFAULT})")
    ap.add_argument("--mip-gap", dest="mip_gap", type=float, default=MIP_GAP_DEFAULT,
                    help=f"Gap relativo del solver (def. {MIP_GAP_DEFAULT})")
    ap.add_argument("--only", nargs="+", default=None, metavar="ESC",
                    help="Resolver sólo estos escenarios: " + " ".join(s for s, _, _, _ in SWEEPS))
    ap.add_argument("--no-lex", dest="no_lex", action="store_true",
                    help="Solve mono-objetivo puro en vez de la cadena lexicográfica "
                         "(más barato; pierde comparabilidad con la referencia)")
    ap.add_argument("--force", action="store_true",
                    help="Re-resolver corridas que ya tienen monoobj_results.json")
    ap.add_argument("--dry-run", dest="dry_run", action="store_true",
                    help="Imprime el plan y los tamaños de modelo, sin resolver")
    args = ap.parse_args()

    bad_obj = [o for o in args.objectives if o not in sens.OBJ_SENSE]
    if bad_obj:
        raise SystemExit(f"Objetivo(s) desconocido(s): {bad_obj}. "
                         f"Opciones: {sorted(sens.OBJ_SENSE)}")
    objectives = list(args.objectives)
    use_lex = not args.no_lex
    REFERENCE_DIR = args.reference

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
        f"results_articulo_{args.zone}_{ts}")

    print(f"\n=== Sensibilidad OFAT sobre soluciones mono-objetivo | zona: {args.zone} ===")
    print(f"  outdir      : {outdir}")
    print(f"  método      : re-optimización "
          + ("con cadena lexicográfica [Mavrotas 2009 §3.1]" if use_lex else "mono-objetivo pura (--no-lex)"))
    print(f"  objetivos   : {objectives}")
    if use_lex:
        for o in objectives:
            print(f"      fila {o}: {' → '.join(sens._lex_row_order(o, PAYOFF_OBJS))}")
    print(f"  solver      : TL={args.time_limit}s  gap={args.mip_gap}")
    n_solves = len(runs) * len(objectives) * (len(PAYOFF_OBJS) if use_lex else 1)
    print(f"  corridas    : {len(runs)} configuraciones × {len(objectives)} objetivos "
          f"≈ {n_solves} solves máx")
    print()

    ref_anchors, ref_params = load_reference(REFERENCE_DIR)
    missing = [o for o in objectives if o not in ref_anchors]
    if missing:
        print(f"  [warning] La referencia no tiene ancla para {missing}: esas filas se "
              f"resolverán en frío (sin warm start ni comparación).")

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
        done = os.path.exists(os.path.join(rundir, "monoobj_results.json"))

        print(f"\n{'='*78}\n[{n}/{len(runs)}] configuración '{run_id}'"
              f"{'  (ya resuelta: se reutiliza)' if done and not args.force else ''}\n{'='*78}")
        diff = "  ".join(f"{k}={v}" for k, v in cfg.items() if BASE_CFG.get(k) != v)
        print(f"  config: {diff if diff else 'caso base (sin desviaciones)'}")

        if done and not args.force:
            status = dict(run_id=run_id, status="reutilizada", error=None, wall_time_s="")
        else:
            status = solve_run(run_id, cfg, args.zone, rundir, ref_anchors,
                               objectives, args.time_limit, args.mip_gap, use_lex)
            print(f"  -> {status['status']}  ({status.get('wall_time_s', '?')}s)")

        loaded = load_run(rundir) or {}
        results[run_id] = dict(status=status,
                               soluciones=loaded.get("soluciones", {}),
                               payoff=loaded.get("payoff", {}),
                               geo=loaded.get("geo", {}))

        write_aggregates(outdir, plan, results, runs, ref_anchors, objectives)

    make_comparison_plots(outdir, plan, results, ref_anchors, objectives)
    write_aggregates(outdir, plan, results, runs, ref_anchors, objectives)

    ok = sum(1 for r in results.values() if r["status"].get("status") in ("ok", "reutilizada"))
    err = [rid for rid, r in results.items() if r["status"].get("status") == "error"]
    print(f"\n{'='*78}")
    print(f"Barrido terminado en {(time.time()-t_all)/60:.1f} min — "
          f"{ok}/{len(runs)} configuraciones con resultados.")
    if err:
        print(f"  con error: {err}  (ver runs/<id>/run.log)")
    print(f"\nTablas agregadas en {outdir}:")
    for f in ("sensibilidad_monoobj.csv", "cardinalidad.csv",
              "plantas_frecuencia.csv", "payoff_tables.csv", "runs_manifest.json"):
        print(f"    {f}")
    print()


if __name__ == "__main__":
    main()
