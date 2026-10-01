#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_instances.py — precómputo de instancias geográficas de Nancy.

Lee el Excel UNA SOLA VEZ y serializa el "núcleo geográfico" (nodos, coords,
arcos de distancia, oferta w_imt y parámetros) de cada zona a JSON, para que
sensibilidad.py no tenga que re-parsear el .xlsx de ~100 MB en cada corrida.

FUENTE ÚNICA: Base_de_datos_FINAL_v2.xlsx (BD consolidada y aumentada). NO se
lee ninguna otra fuente de instancia (ni Base_de_datos_FINAL.xlsx v1 ni
nancy_milp_instance_v2.xlsx): si el archivo v2 no existe, el script aborta.
Cambios de esquema respecto de la BD anterior:
  - Hoja de escalares:  Scalars      → Parameters   (encabezado en fila 0)
  - Polímeros:          Polymers (4) → Polymer (3)  (sin "Otros")
  - Edificios:          Buildings_I  → Buildings     (trae W_1..W_12 mensuales)
  - Oferta de residuos: (calculada)  → w_imt         (i,mes,m ya precomputado)
  - Dist. I→J:          Dist_IJ (km) → D_Build_Hub   (total_cost en metros)
  - Dist. J→K:          Dist_JK (km) → D_Cand_Hub    (facility→hub, metros)
  - Hubs:               Hubs_J       → Hubs          (SIN capacidad Q^J)
  - Plantas:            Facilities_K → Facilities     (class=tipología, con x/y)
  - Estacionalidad:     (hardcodeada)→ ya incluida en W_t / w_imt

Uniones verificadas contra los datos:
  - w_imt."i (fid)" se une por fid a Buildings.fid (conjuntos idénticos).
  - D_Build_Hub.origin_id NO es el fid: es el índice de fila secuencial (1..N)
    del edificio en la hoja Buildings. Se mapea posición→fid por el orden de la
    hoja. destination_id = Hubs.fid (1..82).
  - D_Cand_Hub.origin_id = Facilities.fid (1..46); destination_id = Hubs.fid.
    (Es decir: origen=planta k, destino=hub j → arco J→K = (dest, origin, d).)
  - total_cost está en METROS (validado vs. distancia euclídea) → /1000 = km.

AVISO DE DATOS: la hoja D_Build_Hub se exportó hasta el límite de filas de Excel
(1.048.575 filas de datos). Cubre los edificios secuenciales 1..12.787 completos;
los últimos ~671 edificios quedan sin arcos a hubs y se descartan por
conectividad. Re-exportar D_Build_Hub podado por distancia evitaría la pérdida.

Partición (decisión de tesis): el centro histórico se AÍSLA como un disco
central alrededor del centroide ponderado por residuos; el anillo periférico se
divide en 5 sectores angulares BALANCEADOS. Total = 6 zonas + una 'full'.

Uso:
  python build_instances.py            # reconstruye instances/ desde el Excel
"""

import os, sys, json, math
import numpy as np
import pandas as pd
import openpyxl

try:                                   # consola Windows (cp1252) → UTF-8
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------------------------------------------------------------------------
# FUENTE ÚNICA: Base_de_datos_FINAL_v2.xlsx — BD consolidada y AUMENTADA por
# build_augmented_bd.py con las columnas de capacidad/CAPEX reales por planta
# (Qproc_kg_month, Qpre_kg_month, Ffix_eur), la capacidad real por hub
# (Qhub_kg_month) y el 4º polímero "Other/Mixed" (m=4).
#
# Esta es la ÚNICA fuente de datos de instancia. Se eliminó deliberadamente el
# fallback a Base_de_datos_FINAL.xlsx (v1) y no se lee nancy_milp_instance_v2.xlsx
# ni ninguna otra planilla: si el v2 no está, el script ABORTA en vez de leer una
# fuente distinta en silencio. Todos los parámetros que el código necesita y que
# NO están en el v2 quedan hardcodeados explícitamente abajo (buscar "HARDCODEADO").
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INSTANCE = os.path.join(_SCRIPT_DIR, "..", "data", "raw", "Base_de_datos_FINAL_v2.xlsx")
OUTDIR   = os.path.join(_SCRIPT_DIR, "..", "data", "instances")
if not os.path.exists(INSTANCE):
    raise FileNotFoundError(
        f"No se encontró la fuente única de instancia: {INSTANCE}. "
        "El pipeline lee EXCLUSIVAMENTE Base_de_datos_FINAL_v2.xlsx "
        "(corre build_augmented_bd.py para generarla)."
    )

# Horizonte de arcos a cachear (km). Deben ser ≥ que el mayor D_BAR_* de escenario.
#   Escenarios actuales: D_BAR_IJ ∈ [0.15, 0.30],  D_BAR_JK ∈ [2.0, 3.5].
# HARDCODEADO: no está en Base_de_datos_FINAL_v2 — valor = cota de precómputo de
#   arcos (km), fijada por el mayor radio de elegibilidad D_BAR_* de los escenarios
#   de sensibilidad.py; el v2 sólo trae la distancia por arco, no este horizonte.
ARC_IJ_MAX = 0.50
ARC_JK_MAX = 3.50

# HARDCODEADO: no está en Base_de_datos_FINAL_v2 — valor = decisión de tesis
#   (partición geográfica: centro aislado + 5 sectores angulares balanceados).
N_SECTORS  = 5          # sectores angulares del anillo periférico (+ centro = 6)

# ---------------------------------------------------------------------------
# ESTACIONALIDAD SUAVE (corrección de datos)
# ---------------------------------------------------------------------------
# La BD (Buildings.W_1..W_12 → w_imt) trae un perfil mensual DENTADO e IDÉNTICO
# para los 13.458 edificios:
#     [1.29, 0.43, 0.60, 0.60, 0.60, 1.38, 1.38, 1.38, 1.33, 0.43, 1.29, 1.29]
# con desplomes anómalos de un solo mes en feb/oct (0.43×) y mesetas por bloque
# (mar=abr=may, jun=jul=ago). Eso genera las "V" agudas del margen mensual.
# Aquí se REDISTRIBUYE la masa anual de cada (edificio, polímero) con una curva
# estacional SUAVE (sinusoide con pico en verano), CONSERVANDO EXACTAMENTE el
# total anual por (i, m): solo cambia la FORMA mensual, no la cantidad. Poner
# SMOOTH_SEASONALITY=False reconstruye el perfil dentado original de la BD.
# HARDCODEADO: no está en Base_de_datos_FINAL_v2 — valor = corrección de forma
#   estacional. El v2 SÍ trae la masa mensual (Buildings.W_1..W_12 → w_imt) y una
#   hoja 'Estacionalidad', pero su perfil es dentado (feb/oct=0.43×) e idéntico
#   para todos los edificios; aquí se REEMPLAZA la FORMA mensual por una sinusoide
#   suave conservando EXACTAMENTE el total anual por (i,m) del v2. Estos 3 valores
#   (activar/amplitud/mes-pico) definen esa forma y no provienen de la BD.
SMOOTH_SEASONALITY = True     # False → conserva el perfil dentado original del v2
SEASON_AMPLITUDE   = 0.20     # ±20% alrededor de la media (0 = plano)
SEASON_PEAK_MONTH  = 7        # mes del pico de generación (7 = julio)

def _seasonal_shares(amplitude=SEASON_AMPLITUDE, peak_month=SEASON_PEAK_MONTH):
    """12 pesos mensuales suaves que SUMAN 1 (sinusoide, sin desplomes)."""
    raw = [1.0 + amplitude*math.cos(2.0*math.pi*(t-peak_month)/12.0)
           for t in range(1, 13)]
    s = sum(raw)
    return [r/s for r in raw]

def _smooth_seasonality(wimt):
    """Redistribuye la masa anual de cada (edificio, polímero) con una curva
    suave, conservando el total anual por (i, m). Misma estructura de retorno
    que _read_excel: dict fid → [[m, t, w], ...] con t = 1..12."""
    shares = _seasonal_shares()
    out = {}
    for f, rows in wimt.items():
        annual = {}
        for m, t, w in rows:
            annual[m] = annual.get(m, 0.0) + w
        new_rows = []
        for m, a in annual.items():
            for t in range(1, 13):
                new_rows.append([int(m), int(t), a*shares[t-1]])
        out[f] = new_rows
    return out


# ---------------------------------------------------------------------------
# CAPACIDAD DE HUB (asumida — la BD final no trae Q^J)
# ---------------------------------------------------------------------------
# HARDCODEADO: no está en Base_de_datos_FINAL_v2 — valor = FALLBACK por tipo de
#   establecimiento (Lycée / Collège / École / marché), sólo activo si la hoja Hubs
#   NO trae la columna Qhub_kg_month. El v2 SÍ la trae (capacidad real calculada por
#   build_augmented_bd.py), así que con la fuente actual esta rama queda INACTIVA;
#   se conserva como respaldo para una BD sin esa columna. CALIBRAR antes de publicar.
HUB_CAP_MARKET   = 3000.0   # marché / mercado
HUB_CAP_LYCEE    = 5000.0   # lycée (secundaria alta, gran afluencia)
HUB_CAP_COLLEGE  = 4000.0   # collège (secundaria media)
HUB_CAP_ECOLE    = 2500.0   # école primaire/maternelle/élémentaire
HUB_CAP_DEFAULT  = 3000.0   # tipo desconocido

def _hub_capacity(hub_kind, denom):
    """Capacidad de throughput por hub (kg/período) según tipo de establecimiento."""
    k = str(hub_kind).strip().lower() if hub_kind is not None else ""
    d = str(denom).upper() if denom is not None else ""
    if k == "market":
        return HUB_CAP_MARKET
    if "LYCEE" in d or "LYCÉE" in d:
        return HUB_CAP_LYCEE
    if "COLLEGE" in d or "COLLÈGE" in d or "SEGPA" in d:
        return HUB_CAP_COLLEGE
    if ("ECOLE" in d or "ÉCOLE" in d or "MATERNELLE" in d
            or "PRIMAIRE" in d or "ELEMENTAIRE" in d or "ÉLÉMENTAIRE" in d):
        return HUB_CAP_ECOLE
    return HUB_CAP_DEFAULT


# ---------------------------------------------------------------------------
# DIFERENCIACIÓN POR SITIO (decisión de tesis)
# ---------------------------------------------------------------------------
# HARDCODEADO: no está en Base_de_datos_FINAL_v2 — valor = FALLBACK de capacidad/
#   CAPEX por tipología de sitio (columna `class`), sólo activo si Facilities NO
#   trae Qproc_kg_month/Qpre_kg_month/Ffix_eur. El v2 SÍ las trae (valores reales de
#   build_augmented_bd.py), así que con la fuente actual esta rama queda INACTIVA
#   (ver _k_caps, que prioriza las columnas reales). Se conserva como respaldo para
#   una BD sin esas columnas. CALIBRAR con datos reales por sitio antes de publicar.
SITE_BY_TYPOLOGY = {
    "industrial": dict(Qproc=6000.0, Qpre=3000.0, Ffix=70000.0),
    "works":      dict(Qproc=5000.0, Qpre=2500.0, Ffix=65000.0),
    "brownfield": dict(Qproc=4000.0, Qpre=2000.0, Ffix=95000.0),
    "mixed":      dict(Qproc=3000.0, Qpre=1500.0, Ffix=45000.0),
}
SITE_DEFAULT = dict(Qproc=4000.0, Qpre=2000.0, Ffix=60000.0)

def _site_params(typ):
    """Capacidad/CAPEX por sitio según tipología (cae a SITE_DEFAULT si desconocida)."""
    return SITE_BY_TYPOLOGY.get(str(typ).strip().lower(), SITE_DEFAULT)


def _is_num(v):
    """True si v es un número real (no None, no NaN)."""
    try:
        return v is not None and not math.isnan(float(v))
    except (TypeError, ValueError):
        return False


def _k_caps(r):
    """(Qproc, Qpre, Ffix) de una planta candidata.

    Prioriza los valores REALES de la BD v2 (columnas Qproc/Qpre/Ffix que
    build_augmented_bd.py calcula a partir de la oferta alcanzable por
    distancia). Se usan aunque sean 0.0 (una planta sin oferta alcanzable
    tiene capacidad 0 y, correctamente, el modelo no la abrirá). Solo si las
    columnas NO existen en la BD (valores None → BD antigua) se cae a la
    estimación por tipología (_site_params).
    """
    qproc = getattr(r, "Qproc", None)
    qpre  = getattr(r, "Qpre", None)
    ffix  = getattr(r, "Ffix", None)
    if _is_num(qproc) and _is_num(qpre) and _is_num(ffix):
        return float(qproc), float(qpre), float(ffix)
    sp = _site_params(getattr(r, "typ", None))
    return sp["Qproc"], sp["Qpre"], sp["Ffix"]


# ---------------------------------------------------------------------------
# LECTURA DEL EXCEL (una sola vez, streaming para las hojas grandes)
# ---------------------------------------------------------------------------
def _header_index(ws):
    """Devuelve (idx, header) de la primera fila (encabezado en fila 0)."""
    for row in ws.iter_rows(min_row=1, max_row=1, values_only=True):
        header = list(row)
        return {h: j for j, h in enumerate(header)}, header
    return {}, []

def _all_idx(header):
    """name → [todas las posiciones]. La BD v2 puede traer una columna repetida
    (aumentada dos veces): _header_index sólo conserva la ÚLTIMA, que en este
    archivo quedó en None (se escribió como fórmula sin recalcular). Con todas
    las posiciones podemos tomar el primer valor no nulo (ver _first_val)."""
    d = {}
    for j, h in enumerate(header):
        d.setdefault(h, []).append(j)
    return d

def _first_val(row, idxs):
    """Primer valor no nulo de la fila entre las posiciones dadas (columnas
    duplicadas). Devuelve None si todas están vacías o fuera de rango."""
    for j in idxs:
        if j < len(row) and row[j] is not None:
            return row[j]
    return None

def _to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

def _read_excel(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)

    # --- Parameters (escalares): Symbol → Value ---
    ws = wb["Parameters"]; idx, _ = _header_index(ws)
    sc = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        sym = row[idx["Symbol"]]
        if sym is None:
            continue
        val = _to_float(row[idx["Value"]])
        if val is not None:
            sc[str(sym).strip()] = val

    # --- Polymer (PET/HDPE/PP + Other/Mixed): lectura posicional (encabezados
    # con unicode). Se DEDUPLICA por m: la BD v2 doble-aumentada puede traer la
    # fila m=4 repetida, y una lista con m duplicado haría que sensibilidad.py
    # iterara ese polímero dos veces en las sumas del modelo (doble conteo). ---
    ws = wb["Polymer"]
    pol_rows = []
    _seen_m = set()
    for r_i, row in enumerate(ws.iter_rows(values_only=True)):
        if r_i == 0:
            continue
        if row[0] is None:
            continue
        m_id = int(row[0])
        if m_id in _seen_m:
            continue
        _seen_m.add(m_id)
        pol_rows.append([
            m_id, str(row[1]),
            float(row[2]), float(row[3]), float(row[4]), float(row[5]),
            float(row[6]), float(row[7]), float(row[8]), float(row[9]),
        ])
    pol = pd.DataFrame(pol_rows, columns=["m","poly","zeta","etaJ","etaK",
                                          "price","cpre","cproc","epre","eproc"])

    # --- Buildings: fid, x, y, w_anual = Σ W_1..W_12 (pre-ζ, kg/año) ---
    ws = wb["Buildings"]; idx, _ = _header_index(ws)
    w_cols = [idx[f"W_{k}"] for k in range(1, 13)]
    seq_to_fid = []          # índice secuencial (1-based, orden de hoja) → fid
    res_rows = []
    usage_counts = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        fid = row[idx["fid"]]
        if fid is None:
            continue
        fid = int(fid)
        seq_to_fid.append(fid)                       # todos los edificios, en orden
        x = _to_float(row[idx["x"]]); y = _to_float(row[idx["y"]])
        w_annual = sum((_to_float(row[c]) or 0.0) for c in w_cols)
        u1 = row[idx.get("bdtopo_USAGE1")] if "bdtopo_USAGE1" in idx else None
        usage_counts[u1] = usage_counts.get(u1, 0) + 1
        if x is None or y is None or w_annual <= 0.0:
            continue
        res_rows.append([fid, x, y, w_annual])
    res = pd.DataFrame(res_rows, columns=["i","x","y","w"])
    surviving = set(res["i"].astype(int))
    print(f"  usage (bdtopo_USAGE1): "
          + ", ".join(f"{k}={v}" for k, v in sorted(usage_counts.items(),
                                                     key=lambda kv: -kv[1])[:6]))

    # --- w_imt: oferta precomputada por (fid, mes, polímero) ---
    # Dict fid → {(m, t): w}. Solo edificios que sobreviven (w_anual > 0).
    ws = wb["w_imt"]; idx, _ = _header_index(ws)
    ci, ct, cm, cw = idx["i (fid)"], idx["t (mes)"], idx["m (polímero)"], idx["w_imt"]
    wimt = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        f = row[ci]
        if f is None:
            continue
        f = int(f)
        if f not in surviving:
            continue
        w = _to_float(row[cw])
        if w is None or w <= 0.0:
            continue
        wimt.setdefault(f, []).append([int(row[cm]), int(row[ct]), w])

    # --- Hubs: fid, x, y, Q. Usa Qhub_kg_month real (BD v2, calculado por
    # build_augmented_bd.py a partir de la oferta alcanzable) si la columna
    # existe; si no (BD vieja), cae a la estimación heterogénea por tipo. ---
    ws = wb["Hubs"]; idx, header = _header_index(ws)
    aidx = _all_idx(header)
    has_real_qhub = "Qhub_kg_month" in idx
    hj_rows = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        fid = row[idx["fid"]]
        if fid is None:
            continue
        x = _to_float(row[idx["x"]]); y = _to_float(row[idx["y"]])
        if x is None or y is None:
            continue
        if has_real_qhub:
            Q = _to_float(_first_val(row, aidx["Qhub_kg_month"])) or 0.0
        else:
            kind  = row[idx["hub"]] if "hub" in idx else None
            denom = row[idx["denominatio"]] if "denominatio" in idx else None
            Q = _hub_capacity(kind, denom)
        name = row[idx["appellation"]] if "appellation" in idx else None
        hj_rows.append([int(fid), str(name) if name else f"Hub {int(fid)}",
                        x, y, Q])
    hj = pd.DataFrame(hj_rows, columns=["j","name","x","y","Q"])
    print(f"  Hubs: capacidad {'REAL (Qhub_kg_month, BD v2)' if has_real_qhub else 'asumida por tipo (BD sin Qhub_kg_month)'}")

    # --- Facilities: fid, x, y, tipología (class), nombre, Qproc/Qpre/Ffix.
    # Usa las columnas REALES (Qproc_kg_month/Qpre_kg_month/Ffix_eur, BD v2,
    # calculadas por build_augmented_bd.py a partir de la oferta alcanzable por
    # distancia) si existen; si no (BD vieja), quedan en None y _zone_payload
    # cae a la estimación por tipología (_site_params). ---
    ws = wb["Facilities"]; idx, header = _header_index(ws)
    aidx = _all_idx(header)
    has_real_cap = all(c in idx for c in ("Qproc_kg_month", "Qpre_kg_month", "Ffix_eur"))
    fk_rows = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        fid = row[idx["fid"]]
        if fid is None:
            continue
        x = _to_float(row[idx["x"]]); y = _to_float(row[idx["y"]])
        if x is None or y is None:
            continue
        typ = row[idx["class"]] if "class" in idx else None
        nm  = row[idx["name"]] if "name" in idx else None
        name = str(nm) if nm else f"Site K{int(fid)}"
        if has_real_cap:
            # _first_val: la BD v2 puede traer estas columnas duplicadas (aumentada
            # dos veces); la última quedó en None → se toma la primera con valor.
            qproc = _to_float(_first_val(row, aidx["Qproc_kg_month"])) or 0.0
            qpre  = _to_float(_first_val(row, aidx["Qpre_kg_month"])) or 0.0
            ffix  = _to_float(_first_val(row, aidx["Ffix_eur"])) or 0.0
        else:
            qproc = qpre = ffix = None
        fk_rows.append([int(fid), name, str(typ) if typ else "", x, y, qproc, qpre, ffix])
    fk = pd.DataFrame(fk_rows, columns=["k","name","typ","x","y","Qproc","Qpre","Ffix"])
    print(f"  Facilities: capacidad/CAPEX {'REAL (Qproc/Qpre/Ffix, BD v2)' if has_real_cap else 'asumida por tipología (BD sin columnas nuevas)'}")

    # --- D_Build_Hub: arcos I→J. origin_id=secuencia edificio → fid; dest=hub. ---
    ws = wb["D_Build_Hub"]; idx, _ = _header_index(ws)
    co, cd, cc = idx["origin_id"], idx["destination_id"], idx["total_cost"]
    n_seq = len(seq_to_fid)
    _cmax = max(co, cd, cc)
    dij_rows = []
    lost_trunc = 0
    for row in ws.iter_rows(min_row=2, values_only=True):
        # openpyxl (read_only) recorta las celdas de cola vacías: las filas cuyo
        # total_cost venía vacío (sin ruteo) llegan como tuplas cortas → se saltan.
        if len(row) <= _cmax:
            continue
        o = row[co]
        if o is None:
            continue
        o = int(o)
        if o < 1 or o > n_seq:
            lost_trunc += 1
            continue
        i_fid = seq_to_fid[o - 1]
        if i_fid not in surviving:
            continue
        d = _to_float(row[cc])
        if d is None:
            continue
        d_km = d / 1000.0                              # metros → km
        if d_km > ARC_IJ_MAX:
            continue
        dij_rows.append([i_fid, int(row[cd]), d_km])
    dij = pd.DataFrame(dij_rows, columns=["i","j","d"])

    # --- D_Cand_Hub: arcos J→K. origin=planta k, destino=hub j. ---
    ws = wb["D_Cand_Hub"]; idx, _ = _header_index(ws)
    co, cd, cc = idx["origin_id"], idx["destination_id"], idx["total_cost"]
    _cmax = max(co, cd, cc)
    djk_rows = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        if len(row) <= _cmax:                          # fila ragged (cola vacía recortada)
            continue
        k = row[co]; j = row[cd]
        if k is None or j is None:
            continue
        d = _to_float(row[cc])
        if d is None:
            continue
        d_km = d / 1000.0
        if d_km > ARC_JK_MAX:
            continue
        djk_rows.append([int(j), int(k), d_km])        # (j, k, d)
    djk = pd.DataFrame(djk_rows, columns=["j","k","d"])

    wb.close()
    if lost_trunc:
        print(f"  [aviso] D_Build_Hub: {lost_trunc:,} filas con origin_id fuera de "
              f"rango (truncamiento de Excel) — edificios sin arcos a hubs.")
    return sc, pol, res, wimt, hj, fk, dij, djk


# ---------------------------------------------------------------------------
# PARTICIÓN GEOGRÁFICA (centro + 5 sectores balanceados por residuos)
# ---------------------------------------------------------------------------
_OCTANTS = [("E", 0), ("NE", 45), ("N", 90), ("NO", 135),
            ("O", 180), ("SO", 225), ("S", 270), ("SE", 315)]

def _compass(angle_deg):
    """Rumbo (octante) más cercano a un ángulo dado en grados [0,360)."""
    a = angle_deg % 360.0
    return min(_OCTANTS, key=lambda o: min(abs(a-o[1]), 360-abs(a-o[1])))[0]


def _partition(res):
    """Devuelve (assign_fn, meta) donde assign_fn(x,y)->nombre_zona.

    Reglas geométricas reproducibles:
      - centroide ponderado por w.
      - disco central de radio R_CENTRO con ~1/6 de la masa total de residuos.
      - anillo restante cortado en N_SECTORS cuñas angulares de masa ~igual.
    """
    x = res["x"].to_numpy(); y = res["y"].to_numpy(); w = res["w"].to_numpy()
    cx = float((x*w).sum()/w.sum()); cy = float((y*w).sum()/w.sum())
    r  = np.hypot(x-cx, y-cy)
    w_total = float(w.sum())

    # --- Centro: disco con ~1/6 de la masa (orden por radio creciente) ---
    order = np.argsort(r)
    cum   = np.cumsum(w[order])
    target_c = w_total / 6.0
    n_c = int(np.searchsorted(cum, target_c)) + 1
    n_c = min(n_c, len(order) - 1)
    R_CENTRO = float(r[order][n_c - 1])

    # --- Periferia: cuñas angulares de masa ~igual ---
    peri = r > R_CENTRO
    ang  = (np.degrees(np.arctan2(y - cy, x - cx))) % 360.0
    pa   = ang[peri]; pw = w[peri]
    a_order = np.argsort(pa)
    pa_s, pw_s = pa[a_order], pw[a_order]
    peri_total = float(pw_s.sum())
    cutw = peri_total / N_SECTORS

    # Bordes angulares acumulando masa; el primer borde arranca en el menor ángulo.
    edges = [float(pa_s[0])]
    acc, k = 0.0, 1
    for a_i, w_i in zip(pa_s, pw_s):
        acc += w_i
        if k < N_SECTORS and acc >= k * cutw:
            edges.append(float(a_i)); k += 1
    edges.append(edges[0] + 360.0)   # cierre cíclico

    # Etiquetas únicas por rumbo del centro angular de cada cuña.
    labels, used = [], set()
    for s in range(N_SECTORS):
        mid = (edges[s] + edges[s+1]) / 2.0
        base = _compass(mid); lab = base; n = 2
        while lab in used:
            lab = f"{base}{n}"; n += 1
        used.add(lab); labels.append(lab)

    def assign(px, py):
        if math.hypot(px-cx, py-cy) <= R_CENTRO:
            return "centro"
        a = math.degrees(math.atan2(py-cy, px-cx)) % 360.0
        # localizar la cuña (manejando el wrap del primer borde)
        for s in range(N_SECTORS):
            lo, hi = edges[s], edges[s+1]
            aa = a + 360.0 if a < edges[0] else a
            if lo <= aa < hi:
                return f"sector_{labels[s]}"
        return f"sector_{labels[-1]}"

    meta = dict(centroid=[cx, cy], r_centro=R_CENTRO,
                edges_deg=[e % 360.0 for e in edges[:-1]],
                sector_labels=[f"sector_{l}" for l in labels],
                w_total=w_total)
    zone_names = ["centro"] + [f"sector_{l}" for l in labels]
    return assign, zone_names, meta


# ---------------------------------------------------------------------------
# ENSAMBLADO Y SERIALIZACIÓN POR ZONA
# ---------------------------------------------------------------------------
def _zone_payload(zone, res_z, hj_z, fk_all, dij, djk, sc, pol, wimt):
    """Construye el dict serializable de una zona.

    I y J se particionan por geometría; K se INCLUYE por alcanzabilidad: toda
    planta a la que algún hub de la zona pueda enviar (arco J→K ≤ ARC_JK_MAX).
    Así cada zona es un subproblema MILP autocontenido y factible.
    """
    I_ids = set(res_z["i"].astype(int))
    J_ids = set(hj_z["j"].astype(int))

    ajk = djk[djk["j"].isin(J_ids)].copy()
    K_ids = set(ajk["k"].astype(int))
    fk_z = fk_all[fk_all["k"].isin(K_ids)].copy()
    ajk = ajk[ajk["k"].isin(K_ids)]
    aij = dij[dij["i"].isin(I_ids) & dij["j"].isin(J_ids)].copy()

    # Oferta w_imt de los edificios de la zona: [i, m, t, w]
    w_list = []
    for i in I_ids:
        for (m, t, w) in wimt.get(i, ()):
            w_list.append([int(i), int(m), int(t), float(w)])

    return dict(
        zone=zone,
        meta=dict(n_I=len(I_ids), n_J=len(J_ids), n_K=len(K_ids),
                  w_total=float(res_z["w"].sum()),
                  arc_ij_max=ARC_IJ_MAX, arc_jk_max=ARC_JK_MAX,
                  n_arc_ij=len(aij), n_arc_jk=len(ajk), n_w_imt=len(w_list)),
        scalars={str(k): float(v) for k, v in sc.items()},
        polymers=[[int(r.m), str(r.poly), float(r.zeta), float(r.etaJ),
                   float(r.etaK), float(r.price), float(r.cpre),
                   float(r.cproc), float(r.epre), float(r.eproc)]
                  for r in pol.itertuples()],
        I=[[int(r.i), float(r.x), float(r.y), float(r.w)] for r in res_z.itertuples()],
        J=[[int(r.j), float(r.x), float(r.y), float(r.Q)] for r in hj_z.itertuples()],
        K=[[int(r.k), float(r.x), float(r.y), str(r.name), *_k_caps(r)]
           for r in fk_z.itertuples()],
        A_IJ=[[int(r.i), int(r.j), float(r.d)] for r in aij.itertuples()],
        A_JK=[[int(r.j), int(r.k), float(r.d)] for r in ajk.itertuples()],
        w_imt=w_list,
    )


def build():
    print(f"Leyendo {INSTANCE} (una sola vez)…")
    sc, pol, res, wimt, hj, fk, dij, djk = _read_excel(INSTANCE)
    if SMOOTH_SEASONALITY:
        wimt = _smooth_seasonality(wimt)
        _sh = _seasonal_shares()
        print(f"  [estacionalidad] perfil mensual SUAVIZADO "
              f"(sinusoide amp={SEASON_AMPLITUDE}, pico=mes {SEASON_PEAK_MONTH}); "
              f"masa anual por (i,m) conservada.")
        print(f"     shares×12 = {[round(s*12, 3) for s in _sh]}")
    print(f"  residenciales={len(res)}  hubs={len(hj)}  plantas={len(fk)}  "
          f"polímeros={len(pol)}  arcos_IJ≤{ARC_IJ_MAX}={len(dij)}  "
          f"arcos_JK≤{ARC_JK_MAX}={len(djk)}  w_imt(edif)={len(wimt)}")

    assign, zone_names, meta = _partition(res)
    print(f"  centroide=({meta['centroid'][0]:.0f},{meta['centroid'][1]:.0f})  "
          f"R_centro={meta['r_centro']:.0f} m  cuñas={meta['edges_deg']}")

    res = res.assign(_z=[assign(xx, yy) for xx, yy in zip(res["x"], res["y"])])
    hj  = hj.assign(_z=[assign(xx, yy) for xx, yy in zip(hj["x"], hj["y"])])

    os.makedirs(OUTDIR, exist_ok=True)
    manifest = dict(
        source=os.path.basename(INSTANCE),
        source_assumptions=dict(
            # La oferta viene precomputada en w_imt (post-ζ). El perfil mensual
            # dentado original de la BD se REEMPLAZA por una sinusoide suave
            # (masa anual conservada) cuando SMOOTH_SEASONALITY=True.
            supply_from="w_imt" + (" (estacionalidad suavizada)" if SMOOTH_SEASONALITY else " (estacionalidad original BD)"),
            seasonality_smoothed=bool(SMOOTH_SEASONALITY),
            seasonality_shares_x12=[round(s*12, 4) for s in _seasonal_shares()] if SMOOTH_SEASONALITY else None,
            hub_capacity="heterogénea asumida por tipo (ver _hub_capacity)",
            # Factor de generación de residuos (kg/vivienda/año) usado para
            # generar los datos originales en el Excel. Este valor es crucial
            # para que sensibilidad.py pueda recalcular la oferta con nuevos
            # estimados (ej. W_ESTIMATE=15.0 de la encuesta).
            w_annual_per_dwelling=99.0,
        ),
        partition=meta, zones=[])

    # Zona "full": instancia completa (todos los nodos), también desde caché.
    targets = [("full", res, hj)] + [(z, res[res._z == z], hj[hj._z == z]) for z in zone_names]
    for zone, res_z, hj_z in targets:
        if res_z.empty:
            print(f"  [skip] {zone}: sin zonas residenciales"); continue
        payload = _zone_payload(zone, res_z, hj_z, fk, dij, djk, sc, pol, wimt)
        path = os.path.join(OUTDIR, f"zone_{zone}.json")
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, separators=(",", ":"))
        mm = payload["meta"]
        manifest["zones"].append(dict(name=zone, file=os.path.basename(path), **mm))
        print(f"  → {zone:14s} I={mm['n_I']:5d} J={mm['n_J']:3d} K={mm['n_K']:3d}  "
              f"arcs IJ={mm['n_arc_ij']:5d} JK={mm['n_arc_jk']:4d}  "
              f"w={mm['w_total']/1e6:6.2f}M kg/año  ({os.path.getsize(path)/1e3:.0f} KB)")

    with open(os.path.join(OUTDIR, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1, ensure_ascii=False)
    print(f"\nListo. {len(manifest['zones'])} instancias en {OUTDIR}\\")

    _plot_partition(res, hj, fk, meta, zone_names)


def _plot_partition(res, hj, fk, meta, zone_names):
    """Mapa de las 6 zonas (best-effort; requiere matplotlib)."""
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except Exception:
        print("  (matplotlib no disponible; omito el mapa)"); return

    cx, cy = meta["centroid"]; R = meta["r_centro"]
    cmap = plt.cm.tab10
    zcol = {z: cmap(i % 10) for i, z in enumerate(zone_names)}

    fig, ax = plt.subplots(figsize=(10, 9))
    for z in zone_names:
        sub = res[res._z == z]
        ax.scatter(sub["x"], sub["y"], s=6, color=zcol[z], alpha=0.55, label=z, zorder=2)
    ax.scatter(hj["x"], hj["y"], marker="^", s=70, color="k", zorder=4, label="hubs (J)")
    ax.scatter(fk["x"], fk["y"], marker="s", s=90, facecolor="none",
               edgecolor="red", linewidths=1.4, zorder=5, label="plantas (K)")
    ax.add_patch(Circle((cx, cy), R, fill=False, color="k", ls="--", lw=1.2, zorder=3))
    ax.scatter([cx], [cy], marker="*", s=200, color="gold", edgecolor="k", zorder=6)
    ax.set_aspect("equal"); ax.grid(alpha=0.25)
    ax.set_xlabel("X (L93)"); ax.set_ylabel("Y (L93)")
    ax.set_title("Partición de Nancy: centro + 5 sectores balanceados por residuos")
    ax.legend(loc="upper right", fontsize=8, markerscale=1.5)
    out = os.path.join(OUTDIR, "partition_map.png")
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    print(f"  mapa → {out}")


if __name__ == "__main__":
    build()
