# -*- coding: utf-8 -*-
"""
sensibilidad.py — optimización tri-objetivo de la red de reciclaje de Nancy por
ε-constraint aumentado (AUGMECON, Mavrotas 2009), con el bypass de AUGMECON2
(Mavrotas & Florios 2013) activable por flag.

Configuración pruebafull222 (la que se grilla): Z4=VAN total (primario) con Z2=GEI
y Z3=cobertura como ε-constraints. Z1 (VAN max-min de la peor planta) y Z5 (output
recuperado) quedan COMENTADOS como objetivos —su rango era degenerado/alineado y
colapsaba el frente— pero se siguen calculando y reportando en cada punto como
columnas diagnósticas.

Los tres puntos en que la implementación se aparta del método citado están
declarados y aislados en flags (AUGMECON2_BYPASS, EPS_FEAS_TOL_FRAC,
TIEBREAK_OPEN_W); ver la cabecera "Fidelidad al método citado".

Flujo:
  1. Carga la instancia geográfica precomputada (build_instances.py →
     instances/zone_*.json, generada desde Base_de_datos_FINAL_v2.xlsx). NUNCA se
     re-lee el Excel: para cambiar de BD basta re-correr build_instances.py.
  2. Construye el MILP con apertura dinámica de plantas (y[k,t]) y objetivo Z1
     egalitario (max-min del VAN por planta).
  3. Payoff lexicográfico + grilla ε-constraint secuencial con early-exit.
  4. Guarda TODO el detalle (objetivos + variables de decisión por período) en
     econ_results.json y en CSVs para los mapas territoriales/económicos.

Uso:
  python sensibilidad.py --zone full                 # resuelve y grafica
  python sensibilidad.py --zone full --plot          # solo re-grafica lo guardado
  python sensibilidad.py --list                      # escenarios y zonas
"""

import sys, time, json, os, csv, datetime, warnings
import pandas as pd
import numpy as np
import pyomo.environ as pyo
from scipy.spatial.distance import pdist, squareform
    
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

warnings.filterwarnings("ignore", message=".*Loading a SolverResults.*aborted.*")

# ===========================================================================
# CONFIGURACIÓN
# ===========================================================================
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INSTANCES_DIR = os.path.join(REPO_ROOT, "data", "instances")

ZONE = "full"
HORIZON = 48
ALPHA = 0.00407

# Precio de producto por polímero, p_m [€/kg]. Con USE_PRODUCT_PRICES=True se
# usaban precios de artículo (2.20/4.00/4.00) que daban márgenes de 1.5–2.9 €/kg
# de entrada y hacían el reciclaje INCONDICIONALMENTE rentable → reciclar todo era
# a la vez el óptimo económico, ambiental y de output (objetivos alineados) y el
# frente colapsaba a un punto. Se vuelve a los precios de RECICLADO de la BD
# (PET 0.95, HDPE 0.80, PP 0.55): el margen cae a PET +0.48 / HDPE +0.34 / PP +0.04
# €/kg, PP queda ~en equilibrio y el transporte/costos fijos pasan a decidir →
# reaparece el trade-off económico (VAN vs cobertura).
USE_PRODUCT_PRICES = False
PRICE_PRODUCT = {1: 2.20, 2: 4.00, 3: 4.00}   # override histórico (inflado); inactivo

# Rendimiento de recuperación η_m usado en q = η·x (producto de salida en planta).
# ⚠️ La BD trae DOS columnas que NO deben confundirse:
#   · zeta (0.28/0.32/0.25) = composición del flujo → YA se aplica en el ORIGEN
#     (w_imt = zeta·w_it en build_instances). NO es un rendimiento de planta.
#   · etaK (0.85/0.80/0.75) = rendimiento de recuperación → es el que va en q=η·x.
# Usar zeta como η (comportamiento histórico) aplica la composición DOS veces y
# subestima producto/ingreso/crédito-CO2 ~2.7×, degenerando el frente (max Z1 y
# min Z2 → red vacía (0,0,0)). Por defecto usamos etaK, físicamente correcto.
# El costo unitario c^prod_m se toma directo de la BD (cpre+cproc), SIN MANUF.
YIELD_MODE   = "etaK"    # "etaK" (rendimiento físico, correcto) | "zeta" (composición, antiguo)
YIELD_SCALE  = 1.0       # multiplicador global de η — para el análisis de sensibilidad del yield
ETA_OVERRIDE = None      # dict {m: η} que fuerza valores exactos (tiene prioridad si no es None)

# Crédito de sustitución ε^sub_m [kgCO2eq/kg de salida]. Es el crédito estructural
# del balance neto de GEI (Z2): se resta siempre que el escenario active LCA_CREDIT.
# ⚠️ REDUCIDO a ~50% (LCA NETO tras recolección, sorting, reproceso y pérdidas de
# rendimiento — valor defendible y citable). Con los valores originales (2.5/1.8/2.0)
# Z2/kg quedaba en −1.4…−0.68 (siempre negativo) y el objetivo ambiental seguía
# ALINEADO con el económico. Con estos, Z2/kg CRUZA CERO (PET −0.36, HDPE −0.04,
# PP +0.07): procesar PET reduce CO2 pero procesar PP lo EMITE neto ⇒ el eje
# ambiental por fin se OPONE al VAN y aparece un frente 3D con trade-off real.
# Es el dial ambiental: subir hacia 2.5 vuelve al colapso; bajar refuerza la tensión.
LCA_CREDIT_M = {1: 1.25, 2: 0.90, 3: 1.00, 4: 0.25}
C_MUN = 0.25

# Término de vertedero / contrafactual en Z2 (equivalente al footprint HE^r de
# Ignacio). La masa GENERADA que NO termina como producto q se dispone
# (vertedero/incineración): disposición_m = Γ_m − Σq_m, con factor ε^disp_m
# [kgCO2eq/kg].
# ⚠️ DESACTIVADO: restaba un crédito de vertedero-evitado (2.3–3.1 kgCO2/kg) POR
# ENCIMA del crédito de sustitución LCA por el MISMO kg recuperado → doble conteo
# del beneficio ambiental. Con ambos, Z2 ≈ −3.3 kgCO2/kg (monótono) y el objetivo
# ambiental quedaba pegado al económico. Se deja SOLO el crédito de sustitución
# (LCA_CREDIT_M) como ε^sub único y citable (Z2 ≈ −1.4 kgCO2/kg). Reactivar solo
# si el crédito de sustitución se retira de LCA_CREDIT_M (evitar contarlo dos veces).
USE_DISPOSAL = False
DISPOSAL_EF_M = {1: 2.3, 2: 3.1, 3: 3.1, 4: 1.0}   # factores de disposición; inactivos

# ── Ruptura de colinealidad Z2≡Z3 (min-GEI ≡ max-cobertura) ─────────────────
# El colapso del frente a una recta viene de que min-GEI y max-cobertura son la
# MISMA solución: el crédito neto por kg reciclado (η·(ε^sub+ε^disp) − proceso ≈
# 3.4–3.8 kgCO2/kg) aplasta la emisión de acarreo (eps_tr·1e-3·d ≈ 0.32 kgCO2/kg
# a 1.5 km) ⇒ recolectar SIEMPRE baja Z2. EPS_TR_SCALE escala SÓLO la emisión de
# transporte en Z2 (NO el costo económico ctr), de modo que el material lejano
# deja de ser auto-"verde": con ×15 el cruce de signo cae en ~1.1 km y aparece un
# mínimo interior de Z2 ≠ cobertura máxima. Se fija por escenario (BASE = 1.0,
# físicamente intacto); es el dial a barrer si el frente sigue colineal.
EPS_TR_SCALE = 1.0

SOLVER_TIME_LIMIT = 1000  # corrida NOCTURNA (subido de 1000): hay presupuesto de
                          # toda la noche para cerrar gap por celda. Baja a 180 para
                          # un smoke-test rápido; a 600 para una corrida normal.
SOLVER_MIP_GAP = 0.005  # gap estricto: con 0.03 (~26k€ sobre VAN~887k) el solver
                        # devolvía clones por celda (mismas y[k], mismo punto);
                        # 0.005 separa configuraciones vecinas del frente.
                        # ⚠ La garantía de eficiencia de Mavrotas (2009), Sec. 3.2,
                        # supone solve EXACTO. Con gap>0 las soluciones están dentro
                        # de una banda del óptimo y la garantía no es estricta; la
                        # salvaguarda práctica es el filtro de no-dominancia
                        # post-hoc (_pareto_filter). Este gap también fija la
                        # tolerancia ε de _set_eps (ver EPS_FEAS_TOL_FRAC).
SOLVER_TEE = False
SOLVER_ORDER = ["gurobi", "cbc", "appsi_highs", "glpk"]
# En ε-constraint sólo cambia el RHS (eps2/eps3) entre subproblemas: reusar la
# solución previa como MIP start acelera enormemente y evita los "nosolution"
# (time-limit sin incumbente hallado).
# ⚠ Detalle de implementación del solver, NO una prescripción de Mavrotas (2009):
# el paper no especifica warm-starting. No altera el conjunto de soluciones
# eficientes, sólo el tiempo de cómputo.
SOLVER_WARMSTART = True

RUN_ECONSTRAINT = True
ECON_PRIMARY = "Z4"   # "Z4" = VAN total (Σπ_k). Es el f1 de la ec. (6) de [M09].
# Lista COMPLETA de candidatos a ε-constraint aumentada (todos los objetivos no
# primarios). NO es la lista que necesariamente se grillea: qué objetivo entra a
# la grilla y cuál se fija en su nadir lo decide el RANGO REAL medido en la tabla
# de pagos, no una preferencia editorial. El mecanismo ya existe y no hay que
# reinventarlo — `epsilon_constraint` parte `constrained` en `active` / `loose_only`
# con el umbral (rng[k][1]-rng[k][0]) > 1e-6*max(1,|rng[k][1]|); lo único que hacía
# falta era alimentarlo con una tabla de rangos correcta, y eso es lo que produce
# `lexicographic_payoff_table` (payoff lexicográfico, [M09] Sec. 3.1).
# Un objetivo con rango ~nulo no tiene trade-off que grillar: se fija flojo y se
# sigue calculando/reportando en cada punto (econ_results.json + soluciones.csv).
ECON_CONSTRAINED_ALL = ["Z2", "Z3"]   # config full222: Z1 EXCLUIDO (rango degenerado,
                                      # colapsaba el frente); Z5 (output) queda SÓLO
                                      # como columna diagnóstica (se calcula y reporta
                                      # en cada punto vía `pt`, pero NO se grillea).
ECON_NPTS = 10         # q_i por defecto para los objetivos ausentes de ECON_NPTS_MAP.
# Densidad de grilla POR OBJETIVO (q_i), no uniforme. [M09] Sec. 4: "we can control
# the density of the efficient set representation ... by properly assigning the
# values to the q_i ... trade off between the density of the efficient set and the
# computation time". Con hasta 4 objetivos activos el coste Π q_i crece rápido; la
# respuesta del método citado es bajar la resolución de los ejes menos informativos,
# NO excluirlos de la grilla. El último de `active` es el bucle interno (early-exit
# + bypass), así que conviene darle el q_i más alto.
ECON_NPTS_MAP = {"Z2": 15, "Z3": 15}   # grilla densa pruebafull222: 15×15 nominal
                                        # (225). early-exit + bypass + warmstart la
                                        # reducen a ~40-60 solves reales. pruebafull222
                                        # corrió 12×12; esto es la versión más densa.
ECON_AUGW = 1e-3      # eps de la ec. (6) de Mavrotas (2009); el paper recomienda
                      # un valor "usually between 1e-3 and 1e-6".

# ── Fidelidad al método citado ──────────────────────────────────────────────
# Referencias:
#   [M09] Mavrotas, G. (2009). "Effective implementation of the ε-constraint
#         method in Multi-Objective Mathematical Programming problems".
#         Applied Mathematics and Computation 213:455-465.  Ecs. (1)-(9), Secs. 2-4.
#   [MF13] Mavrotas, G., Florios, K. (2013). "An improved version of the augmented
#         ε-constraint method (AUGMECON2) for finding the exact Pareto set in
#         Multi-Objective Integer Programming problems".
#         Applied Mathematics and Computation 219:9652-9669.
#
# Los tres flags de abajo son los ÚNICOS puntos donde la implementación se aparta
# de [M09] tal cual está escrito. Cada uno se declara, se cita y se puede apagar.

# (a) BYPASS por magnitud de holgura. NO está en [M09]: la Sec. 3.3 / Fig. 4 de
#     [M09] sólo describe early-exit por INFACTIBILIDAD. Saltar puntos de grilla
#     FACTIBLES usando el tamaño de la holgura del objetivo interno es la mejora
#     de [MF13] (coeficiente de bypass b = floor(s_in / step_in)).
#     True  -> AUGMECON2 [MF13] (más rápido, mismo conjunto de Pareto).
#     False -> AUGMECON puro [M09] (sólo early-exit por infactibilidad).
#     Si se deja en True, el manuscrito DEBE citar [MF13], no sólo [M09].
AUGMECON2_BYPASS = True

# (b) TOLERANCIA DE FACTIBILIDAD sobre el RHS ε. NO está en [M09]: la ec. (3) fija
#     el RHS exacto. Es una necesidad numérica de trabajar con MILP (y no LP): las
#     celdas de la tabla de pagos sólo están certificadas dentro del gap relativo
#     del solver, así que fijar eps_k EXACTAMENTE en el óptimo individual de k
#     puede volver infactible el extremo ideal del frente. Se relaja eps_k en
#         ftol_k = EPS_FEAS_TOL_FRAC * SOLVER_MIP_GAP * R_k
#     (R_k = amplitud del rango de k en el payoff = escala natural de la grilla).
#     Queda ATADA al gap real del solver en vez de a una constante suelta. Con los
#     valores por defecto (1e-3 × 0.005) da ftol = 5e-6·R_k, tres órdenes de
#     magnitud por debajo del paso de grilla R_k/(npts-1) ⇒ absorbe el error de
#     optimalidad del payoff sin correr ningún punto de la grilla a otra celda.
EPS_FEAS_TOL_FRAC = 1e-3

# (c) DESEMPATE por nº de plantas abiertas. NO es parte de AUGMECON: la ec. (6) de
#     [M09] es max(f1 + eps*Σ s_k/r_k), sin tercer término. Un peso > 0 aquí
#     penaliza abrir plantas y sesga TODO el frente hacia redes más pequeñas.
#     Se deja en 0.0 (objetivo = ec. (6) textual). Sólo subirlo si hace falta
#     romper empates entre configuraciones con idéntico valor de objetivo y
#     distinto nº de plantas — y en ese caso REPORTARLO como criterio de
#     desempate, explícitamente fuera del método citado.
TIEBREAK_OPEN_W = 0.0

MAKE_PLOTS = True
EXPORT_ZONE_DETAIL = True

# ── Auditoría de gap / curva de optimalidad ─────────────────────────────────
# Tolerancia de gap (%) por encima de la cual una celda se considera "problemática"
# (no cerró): entra en la lista de entrada del modo --time-curve. Las celdas en
# 'timelimit' también se marcan aunque su gap sea bajo.
GAP_AUDIT_TOL_PCT = 100.0
# Modo --time-curve: por cada celda problemática se resuelve secuencialmente con
# estos límites de tiempo crecientes (s), reusando el incumbente previo como
# warmstart. Corte temprano si el gap baja de TIME_CURVE_GAP_STOP_PCT.
TIME_CURVE_LIMITS = [1500,2500,3600]
TIME_CURVE_GAP_STOP_PCT = 0.5

# ===========================================================================
# ESCENARIOS
# ===========================================================================
SCENARIO = "BASE"

_SCENARIOS = {
    # Caso base del manuscrito (§Instance / Apéndice A): radios de elegibilidad
    # d_bar^IJ=0.25 km (250 m) y d_bar^JK=1.5 km, espaciamiento mínimo
    # d_bar_kk=0.6 km, crédito de sustitución estructural activo, δ=θ=β_min=0.
    # THETA reactivado al valor de la BD (500). ⚠️ Con Z3 = output_total [kg] − THETA·Σy,
    # y output ~10^6 kg vs 500·(≈16 plantas)=8k, THETA=500 es ~2 órdenes menor que el
    # output: apenas decorrela Z3. El problema de fondo es que Z3 mide OUTPUT (alineado
    # con VAN), no un objetivo social; decorrelarlo de verdad requiere replantear Z3.
    "BASE": dict(tag="base", D_BAR_IJ=0.25, D_BAR_JK=1.5, D_BAR_KK_KM=0.6,
                 D_BAR_CENTER_K=None,
                 N_CAP=None, THETA=500.0, PRICE_FACTOR=1.0, WASTE_FACTOR=1.0,
                 LCA_CREDIT=True, BETA_MIN=0.0, MUNICIPAL_CREDIT_SHARE=0.25),
    "LCA": dict(tag="lca", D_BAR_IJ=0.30, D_BAR_JK=3.5, D_BAR_KK_KM=0.3,
                D_BAR_CENTER_K=None,
                N_CAP=None, THETA=0.0, PRICE_FACTOR=1.0, WASTE_FACTOR=1.0,
                LCA_CREDIT=True, BETA_MIN=0.0, MUNICIPAL_CREDIT_SHARE=0.0),
    "WIDE": dict(tag="wide", D_BAR_IJ=0.50, D_BAR_JK=3.5, D_BAR_KK_KM=0.3,
                 D_BAR_CENTER_K=None,
                 N_CAP=None, THETA=0.0, PRICE_FACTOR=1.0, WASTE_FACTOR=1.0,
                 LCA_CREDIT=True, BETA_MIN=0.0, MUNICIPAL_CREDIT_SHARE=0.0),
    # Ruptura de colinealidad por ACARREO (palanca elegida): EPS_TR_SCALE sube la
    # emisión de transporte en Z2 para que el material lejano no sea auto-"verde"
    # (min-GEI deja de ser cobertura máxima), y BETA_MIN exige una cobertura piso
    # que mata el extremo económico "flaco" del 55% ⇒ empuja puntos al interior.
    # Créditos LCA/vertedero INTACTOS; radios = BASE. Diales para empujar más:
    # subir EPS_TR_SCALE (20–30) y/o bajar D_BAR_JK a 1.2–1.3 (riesgo: desconectar
    # zonas lejanas, ver el raise de "Instancia desconectada").
    "DECOUP": dict(tag="decoup", D_BAR_IJ=0.25, D_BAR_JK=1.5, D_BAR_KK_KM=0.6,
                   D_BAR_CENTER_K=None,
                   N_CAP=None, THETA=0.0, PRICE_FACTOR=1.0, WASTE_FACTOR=1.0,
                   LCA_CREDIT=True, BETA_MIN=0.3, MUNICIPAL_CREDIT_SHARE=0.0,
                   EPS_TR_SCALE=15.0),
}

if SCENARIO not in _SCENARIOS:
    raise ValueError(f"SCENARIO '{SCENARIO}' no reconocido. Opciones: {list(_SCENARIOS)}")

_S = _SCENARIOS[SCENARIO]
D_BAR_IJ = _S["D_BAR_IJ"]
D_BAR_JK = _S["D_BAR_JK"]
D_BAR_KK_KM = _S.get("D_BAR_KK_KM")
D_BAR_CENTER_K = _S.get("D_BAR_CENTER_K")
N_CAP = _S["N_CAP"]
THETA = _S["THETA"]
PRICE_FACTOR = _S["PRICE_FACTOR"]
WASTE_FACTOR = _S["WASTE_FACTOR"]
LCA_CREDIT = _S["LCA_CREDIT"]
BETA_MIN = _S.get("BETA_MIN", 0.0)
MUNICIPAL_CREDIT_SHARE = _S.get("MUNICIPAL_CREDIT_SHARE", 0.0)
EPS_TR_SCALE = _S.get("EPS_TR_SCALE", 1.0)   # multiplicador de la emisión de acarreo (Z2)

_h_eff = HORIZON if HORIZON is not None else 60
_MONTHS = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]
MONTH_NAMES = (_MONTHS[:_h_eff] if _h_eff <= 12
               else [f"{_MONTHS[t % 12]}A{t // 12 + 1}" for t in range(_h_eff)])

_ts = datetime.datetime.now().strftime("%Y%m%d_%H%M")
ECON_OUTDIR = os.path.join(REPO_ROOT, "results",
                           f"results_sens_{ZONE}_{_S['tag']}_{_ts}")


# ===========================================================================
# CARGA DE DATOS
# ===========================================================================
def _filter_k_by_centrality(res_df, fk_df, max_dist_km):
    if max_dist_km is None or max_dist_km <= 0:
        return fk_df

    if res_df.empty:
        print("  [warning] No hay zonas residenciales para calcular el centroide. No se filtra por centralidad.")
        return fk_df

    # Calcular centroide ponderado por residuos de las zonas residenciales (I)
    weights = res_df["w"].replace(0, 1e-9)
    total_weight = weights.sum()
    if total_weight == 0:
        print("  [warning] Peso total de residuos es cero. No se filtra por centralidad.")
        return fk_df

    cx = (res_df["x"] * weights).sum() / total_weight
    cy = (res_df["y"] * weights).sum() / total_weight

    # Calcular distancia de cada planta k al centroide (coordenadas en metros -> km)
    distances_km = np.hypot(fk_df["x"] - cx, fk_df["y"] - cy) / 1000.0

    # Filtrar plantas
    original_k_count = len(fk_df)
    fk_df_filtered = fk_df[distances_km <= max_dist_km].copy()
    new_k_count = len(fk_df_filtered)

    if new_k_count < original_k_count:
        print(f"  [info] Filtrado por centralidad (D_BAR_CENTER_K={max_dist_km} km): "
              f"{original_k_count} -> {new_k_count} plantas candidatas.")

    return fk_df_filtered


def _load_raw_data(zone):
    path = os.path.join(INSTANCES_DIR, f"zone_{zone}.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"No existe {path}.\nGenera la caché primero: python build_instances.py")
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)

    sc = {k: float(v) for k, v in data["scalars"].items()}
    pol = pd.DataFrame(data["polymers"], columns=["m", "poly", "zeta", "etaJ", "etaK", "price", "cpre", "cproc", "epre", "eproc"])
    res_df = pd.DataFrame(data["I"], columns=["i", "x", "y", "w"])
    hj_df = pd.DataFrame(data["J"], columns=["j", "x", "y", "Q"])
    fk_df = pd.DataFrame(data["K"], columns=["k", "x", "y", "name", "Qproc", "Qpre", "Ffix"])
    dij_full_df = pd.DataFrame(data["A_IJ"], columns=["i", "j", "d"])
    djk_full_df = pd.DataFrame(data["A_JK"], columns=["j", "k", "d"])

    w_imt = {}
    for i, m, t, w in data.get("w_imt", []):
        w_imt.setdefault(int(i), {})[(int(m), int(t))] = float(w)

    if res_df.empty:
        print("  [warning] Zona sin zonas residenciales.")
        return None
    return sc, pol, res_df, hj_df, fk_df, dij_full_df, djk_full_df, w_imt


def _filter_connected_graph(res_df, hj_df, fk_df, dij_full_df, djk_full_df):
    i_nodes, j_nodes, k_nodes = set(res_df.i), set(hj_df.j), set(fk_df.k)

    ajk_df = djk_full_df[(djk_full_df.j.isin(j_nodes)) & (djk_full_df.k.isin(k_nodes)) & (djk_full_df.d <= D_BAR_JK)]
    K = sorted(set(ajk_df.k))
    j_connected_to_k = set(ajk_df.j)

    aij_df = dij_full_df[(dij_full_df.i.isin(i_nodes)) & (dij_full_df.j.isin(j_connected_to_k)) & (dij_full_df.d <= D_BAR_IJ)]
    J = sorted(set(aij_df.j))
    I = sorted(set(aij_df.i))

    if not I or not J or not K:
        raise ValueError(f"Instancia desconectada con D_BAR_IJ={D_BAR_IJ}, D_BAR_JK={D_BAR_JK}.")

    ajk_df = ajk_df[ajk_df.j.isin(J)]
    K = sorted(set(ajk_df.k))
    if not K:
        raise ValueError("El grafo se desconectó al filtrar J: no quedan plantas (K).")

    res_df = res_df[res_df.i.isin(I)].copy()
    hj_df = hj_df[hj_df.j.isin(J)].copy()
    fk_df = fk_df[fk_df.k.isin(K)].copy()
    return I, J, K, res_df, hj_df, fk_df, aij_df, ajk_df


def _build_final_model_dict(I, J, K, sc, pol, res_df, hj_df, fk_df, aij_df, ajk_df, w_imt):
    Q_prod = {int(r.k): float(r.Qproc) for r in fk_df.itertuples()}
    F_fix = {int(r.k): float(r.Ffix) for r in fk_df.itertuples()}

    lam = sc["lambda"] * WASTE_FACTOR
    ctr = sc["c_tr"]
    eptr = sc["eps_tr"] * EPS_TR_SCALE   # escala SÓLO la emisión de acarreo (Z2), no el costo ctr
    alpha = ALPHA if ALPHA is not None else sc.get("alpha", sc.get("alpha_m"))
    horizon = HORIZON if HORIZON is not None else int(sc.get("n_periods", 60))

    A_IJ = {(int(r.i), int(r.j)): float(r.d) for r in aij_df.itertuples()}
    A_JK = {(int(r.j), int(r.k)): float(r.d) for r in ajk_df.itertuples()}

    M_idx, name_m, eta = [], {}, {}
    price, cprod, eprod = {}, {}, {}
    for r in pol.itertuples():
        mid = int(r.m)
        M_idx.append(mid)
        name_m[mid] = r.poly
        if ETA_OVERRIDE is not None:
            base_eta = float(ETA_OVERRIDE.get(mid, r.etaK))
        elif YIELD_MODE == "zeta":
            base_eta = float(r.zeta)
        else:                                    # "etaK": rendimiento físico (correcto)
            base_eta = float(r.etaK)
        eta[mid] = YIELD_SCALE * base_eta
        base_price = PRICE_PRODUCT.get(mid, float(r.price)) if USE_PRODUCT_PRICES else float(r.price)
        price[mid] = base_price * PRICE_FACTOR
        cprod[mid] = float(r.cpre) + float(r.cproc)          # c^prod_m = pre+proc (BD), sin MANUF
        eprod[mid] = float(r.epre) + float(r.eproc)

    TT = list(range(1, horizon + 1))
    w_supply = {}
    for i in I:
        wi = w_imt.get(i, {})
        for mm in M_idx:
            for t in TT:
                month_t = (t - 1) % 12 + 1
                w_supply[(i, mm, t)] = WASTE_FACTOR * wi.get((mm, month_t), 0.0)

    Gamma = {i: sum(w_supply.get((i, mm, t), 0.0) for mm in M_idx for t in TT) for i in I}
    I = sorted(i for i in I if Gamma.get(i, 0.0) > 1e-6)
    if not I:
        raise ValueError("No hay zonas de generación válidas.")
    Gamma = {i: Gamma[i] for i in I}
    w_supply = {(i, mm, t): w for (i, mm, t), w in w_supply.items() if i in I}
    A_IJ = {(i, j): d for (i, j), d in A_IJ.items() if i in I}

    P = {i: Gamma[i] / lam if lam > 0 else 0.0 for i in I}
    Qhub = {int(r.j): float(r.Q) for r in hj_df.itertuples()}
    disc = {t: (1.0 + alpha) ** (-t) for t in TT}

    adj_ij = {i: [] for i in I}
    adj_ji = {j: [] for j in J}
    for (i, j) in A_IJ:
        if i in adj_ij: adj_ij[i].append(j)
        if j in adj_ji: adj_ji[j].append(i)
    adj_jk = {j: [] for j in J}
    adj_kj = {k: [] for k in K}
    for (j, k) in A_JK:
        if j in adj_jk: adj_jk[j].append(k)
        if k in adj_kj: adj_kj[k].append(j)

    disc_sum = sum(disc.values())
    max_margin = max((price[mm] * eta[mm] + C_MUN * MUNICIPAL_CREDIT_SHARE) for mm in M_idx) if M_idx else 0.0
    max_Qprod = max(Q_prod.values()) if Q_prod else 0.0
    M_pi_global = 1.05 * disc_sum * max_Qprod * max_margin
    if M_pi_global <= 0:
        M_pi_global = 1e6

    coord_i = {int(r.i): (float(r.x), float(r.y), float(r.w)) for r in res_df.itertuples() if int(r.i) in I}
    coord_j = {int(r.j): (float(r.x), float(r.y)) for r in hj_df.itertuples() if int(r.j) in J}
    coord_k = {int(r.k): (float(r.x), float(r.y), str(r.name)) for r in fk_df.itertuples() if int(r.k) in K}

    return dict(
        I=I, J=J, K=K, M=M_idx, TT=TT, name=name_m,
        A_IJ=A_IJ, A_JK=A_JK,
        adj_ij=adj_ij, adj_ji=adj_ji, adj_jk=adj_jk, adj_kj=adj_kj,
        w=w_supply, P=P, Gamma=Gamma,
        eta=eta, price=price, cprod=cprod, eprod=eprod,
        ctr=ctr, eptr=eptr, Qhub=Qhub,
        Q_prod=Q_prod, F_fix=F_fix,
        disc=disc, M_pi_global=M_pi_global,
        N_cap=(len(K) if N_CAP is None else N_CAP),
        beta_min=BETA_MIN,
        coord_i=coord_i, coord_j=coord_j, coord_k=coord_k,
        horizon=horizon, alpha=alpha,
        municipal_credit_share=MUNICIPAL_CREDIT_SHARE,
    )


def _cluster_k_nodes(fk_df, d_kk_km_threshold):
    if d_kk_km_threshold is None or d_kk_km_threshold <= 0:
        return fk_df.copy(), {int(k): int(k) for k in fk_df["k"]}

    print(f"  [info] Agrupando plantas candidatas con D_BAR_KK_KM = {d_kk_km_threshold} km...")
    nodes = fk_df.copy()
    d_kk_m = d_kk_km_threshold * 1000.0

    parent = {k: k for k in nodes["k"]}
    def find(k):
        if parent[k] == k:
            return k
        parent[k] = find(parent[k])
        return parent[k]
    def union(k1, k2):
        r1, r2 = find(k1), find(k2)
        if r1 != r2:
            if r1 < r2: parent[r2] = r1
            else: parent[r1] = r2

    coords = nodes[["k", "x", "y"]].set_index("k")
    ks = coords.index.values
    dist_matrix = squareform(pdist(coords.values))
    pairs = [(dist_matrix[a, b], ks[a], ks[b]) for a in range(len(ks)) for b in range(a + 1, len(ks))]
    pairs.sort()
    for dist, k1, k2 in pairs:
        if dist > d_kk_m:
            break
        union(k1, k2)

    final_clusters = {}
    for k in nodes["k"]:
        final_clusters.setdefault(find(k), []).append(k)
    if len(final_clusters) == len(nodes):
        print("  [info] No se realizaron agrupaciones.")
        return fk_df.copy(), {int(k): int(k) for k in fk_df["k"]}

    nodes.set_index("k", inplace=True)
    clustered_rows = []
    for root_k, members in final_clusters.items():
        cdf = nodes.loc[members]
        weights = cdf["Qproc"].replace(0, 1e-6)
        total_w = weights.sum()
        clustered_rows.append({
            "k": root_k,
            "x": (cdf["x"] * weights).sum() / total_w,
            "y": (cdf["y"] * weights).sum() / total_w,
            "Qproc": cdf["Qproc"].sum(),
            "Qpre": cdf["Qpre"].sum(),
            "Ffix": cdf["Ffix"].sum(),
            "name": f"Cluster k={root_k} (+{len(members)-1})" if len(members) > 1 else cdf.iloc[0]["name"],
        })
    new_fk_df = pd.DataFrame(clustered_rows).sort_values("k").reset_index(drop=True)
    print(f"  [info] {len(fk_df)} plantas reducidas a {len(new_fk_df)} por agrupación.")
    return new_fk_df, {k: find(k) for k in nodes.index}


def load_instance(zone=None):
    zone = zone or ZONE
    raw_data = _load_raw_data(zone)
    if raw_data is None:
        return None
    sc, pol, res_df, hj_df, fk_df, dij_full_df, djk_full_df, w_imt = raw_data

    # Filtrar plantas candidatas por distancia al centroide de generación
    if D_BAR_CENTER_K:
        fk_df = _filter_k_by_centrality(res_df, fk_df, D_BAR_CENTER_K)

    fk_df_clustered, aggregation_map = _cluster_k_nodes(fk_df, D_BAR_KK_KM)
    if len(aggregation_map) > len(fk_df_clustered):
        djk_full_df["k"] = djk_full_df["k"].map(aggregation_map)
        djk_full_df = djk_full_df.groupby(["j", "k"], as_index=False).agg({"d": "min"})
        fk_df = fk_df_clustered

    I, J, K, res_df, hj_df, fk_df, aij_df, ajk_df = _filter_connected_graph(
        res_df, hj_df, fk_df, dij_full_df, djk_full_df)
    return _build_final_model_dict(I, J, K, sc, pol, res_df, hj_df, fk_df, aij_df, ajk_df, w_imt)


# ===========================================================================
# MODELO
# ===========================================================================
def build_model(D):
    m = pyo.ConcreteModel("nancy_3obj_dynamic_egalitarian")
    I, J, K, MM, TT = D["I"], D["J"], D["K"], D["M"], D["TT"]
    AIJ = list(D["A_IJ"].keys())
    AJK = list(D["A_JK"].keys())

    eta, price, cprod, eprod = D["eta"], D["price"], D["cprod"], D["eprod"]
    ctr, eptr = D["ctr"], D["eptr"]
    dij, djk = D["A_IJ"], D["A_JK"]
    aji, ajk_d, akj = D["adj_ji"], D["adj_jk"], D["adj_kj"]
    aij_d = D["adj_ij"]
    disc, w = D["disc"], D["w"]
    T_last = TT[-1]

    m.y = pyo.Var(K, TT, domain=pyo.Binary)
    m.delta = pyo.Var(K, TT, domain=pyo.NonNegativeReals, bounds=(0, 1))
    m.b = pyo.Var(AIJ, MM, TT, domain=pyo.NonNegativeReals)
    m.x = pyo.Var(AJK, MM, TT, domain=pyo.NonNegativeReals)
    m.q = pyo.Var(K, MM, TT, domain=pyo.NonNegativeReals)
    m.phi = pyo.Var(domain=pyo.Reals)

    _prev = {TT[i]: TT[i - 1] for i in range(1, len(TT))}

    def _delta_rule(m, k, t):
        return m.delta[k, t] == (m.y[k, t] if t == TT[0] else m.y[k, t] - m.y[k, _prev[t]])
    m.delta_def = pyo.Constraint(K, TT, rule=_delta_rule)

    def _mono_rule(m, k, t):
        if t == TT[0]:
            return pyo.Constraint.Skip
        return m.y[k, t] >= m.y[k, _prev[t]]
    m.monotone = pyo.Constraint(K, TT, rule=_mono_rule)

    def _one_open_event_rule(m, k):
        return sum(m.delta[k, t] for t in TT) <= 1
    m.one_open_event = pyo.Constraint(K, rule=_one_open_event_rule)

    def _zone_supply_rule(m, i, mm, t):
        if not aij_d.get(i):
            return pyo.Constraint.Skip
        return sum(m.b[i, j, mm, t] for j in aij_d[i]) <= w[i, mm, t]
    m.zone_supply = pyo.Constraint(I, MM, TT, rule=_zone_supply_rule)

    def _hub_cap_rule(m, j, t):
        if not D["adj_ji"].get(j):
            return pyo.Constraint.Skip
        return sum(m.b[i, j, mm, t] for i in D["adj_ji"][j] for mm in MM) <= D["Qhub"][j]
    m.hub_cap = pyo.Constraint(J, TT, rule=_hub_cap_rule)

    def _hub_balance_rule(m, j, mm, t):
        return sum(m.b[i, j, mm, t] for i in aji[j]) == sum(m.x[j, k, mm, t] for k in ajk_d[j])
    m.hub_balance = pyo.Constraint(J, MM, TT, rule=_hub_balance_rule)

    m.hub_active = pyo.Constraint(J, TT,
        rule=lambda m, j, t: sum(m.b[i, j, mm, t] for i in aji[j] for mm in MM)
        <= D["Qhub"][j] * sum(m.y[k, t] for k in ajk_d[j]))

    m.fac_balance = pyo.Constraint(K, MM, TT,
        rule=lambda m, k, mm, t: m.q[k, mm, t] == eta[mm] * sum(m.x[j, k, mm, t] for j in akj[k]))

    m.intake_cap = pyo.Constraint(K, TT,
        rule=lambda m, k, t: sum(m.x[j, k, mm, t] for j in akj[k] for mm in MM) <= D["Q_prod"][k] * m.y[k, t])

    m.cardinality = pyo.Constraint(expr=sum(m.y[k, T_last] for k in K) <= D["N_cap"])

    if D["beta_min"] > 0:
        _all_ij = list(D["A_IJ"].keys())
        _gen_t = {t: sum(D["w"][i, mm, t] for i in I for mm in MM) for t in TT}
        m.beta_min_c = pyo.Constraint(TT,
            rule=lambda m, t: sum(m.b[i, j, mm, t] for (i, j) in _all_ij for mm in MM) >= D["beta_min"] * _gen_t[t])

    _arcs_jk = [(j, k) for j, ks in D["adj_jk"].items() for k in ks]

    def _pi_rule(m, k):
        revenue = sum(disc[t] * sum(price[mm] * m.q[k, mm, t] for mm in MM) for t in TT)
        mun_credit = C_MUN * D["municipal_credit_share"] * sum(
            disc[t] * sum(m.x[j, k, mm, t] for j in akj[k] for mm in MM) for t in TT)
        transport = ctr * 1e-3 * sum(
            disc[t] * sum(djk[(j, k)] * m.x[j, k, mm, t] for j in akj[k] for mm in MM) for t in TT)
        processing = sum(
            disc[t] * sum(cprod[mm] * m.x[j, k, mm, t] for j in akj[k] for mm in MM) for t in TT)
        fixed_cost = sum(D["F_fix"][k] * disc[t] * m.delta[k, t] for t in TT)
        return revenue + mun_credit - transport - processing - fixed_cost
    m.pi = pyo.Expression(K, rule=_pi_rule)

    m.phi_constraint = pyo.Constraint(K, rule=lambda m, k: m.phi <= m.pi[k] + D["M_pi_global"] * (1 - m.y[k, T_last]))
    m.Z1 = pyo.Expression(expr=m.phi)
    m.phi_if_closed = pyo.Constraint(expr=m.phi <= D["M_pi_global"] * sum(m.y[k, T_last] for k in K))

    m.NPV = pyo.Expression(rule=lambda m: sum(m.pi[k] for k in K))

    def _Z2_rule(m):
        E = sum(
            eptr * 1e-3 * sum(djk[(j, k)] * m.x[j, k, mm, t] for (j, k) in _arcs_jk for mm in MM)
            + sum(eprod[mm] * m.x[j, k, mm, t] for (j, k) in _arcs_jk for mm in MM)
            for t in TT)
        if LCA_CREDIT:
            E -= sum(LCA_CREDIT_M.get(mm, 0.0) * m.q[k, mm, t] for k in K for mm in MM for t in TT)
        if USE_DISPOSAL:
            # Disposición de lo no recuperado = Γ_m − Σq_m. El término Γ_m es
            # constante (emisión de "no hacer nada"); −ε^disp·Σq premia recuperar.
            gamma_m = {}
            for (ii, mm2, tt2), wv in D["w"].items():
                gamma_m[mm2] = gamma_m.get(mm2, 0.0) + wv
            E += sum(DISPOSAL_EF_M.get(mm, 0.0) * gamma_m.get(mm, 0.0) for mm in MM)
            E -= sum(DISPOSAL_EF_M.get(mm, 0.0) * m.q[k, mm, t] for k in K for mm in MM for t in TT)
        return E
    m.Z2 = pyo.Expression(rule=_Z2_rule)

    m.coverage = pyo.Expression(expr=sum(
        (D["P"][i] / max(D["Gamma"][i], 1e-6)) * sum(m.b[i, j, mm, t] for j in aij_d[i] for mm in MM for t in TT)
        for i in I))
    m.Z3 = pyo.Expression(expr=m.coverage)  # Objetivo 3: cobertura social (paper)

    # Objetivo 5: output total recuperado [kg], a MAXIMIZAR.
    m.total_output = pyo.Expression(expr=sum(m.q[k, mm, t] for k in K for mm in MM for t in TT))
    m.Z5 = pyo.Expression(expr=m.total_output)

    m._D = D
    return m


# ===========================================================================
# SOLVER
# ===========================================================================
def _configure_solver(s, name):
    tl, gap = float(SOLVER_TIME_LIMIT), float(SOLVER_MIP_GAP)
    if name == "gurobi":
        s.options["TimeLimit"] = tl
        s.options["MIPGap"] = gap
        s.options["Presolve"] = 2
        s.options["Cuts"] = 2
        s.options["MIPFocus"] = 1
        s.options["Heuristics"] = 0.2
        s.options["Threads"] = 0
        s.options["DualReductions"] = 0
    elif name == "cbc":
        s.options["seconds"] = tl
        s.options["ratioGap"] = gap
    elif name == "appsi_highs":
        s.options["time_limit"] = tl
        s.options["mip_rel_gap"] = gap
    elif name == "glpk":
        s.options["tmlim"] = int(tl)
        s.options["mipgap"] = gap
    return s


_SOLVER_NAME = None
_WARMSTART_OK = {"gurobi", "cplex", "cbc"}   # backends con MIP-start vía archivo


def make_solver():
    global _SOLVER_NAME
    for name in SOLVER_ORDER:
        try:
            s = pyo.SolverFactory(name)
            if s is not None and s.available(exception_flag=False):
                _SOLVER_NAME = name
                print(f"  solver: {name}  (warmstart={'sí' if SOLVER_WARMSTART and name in _WARMSTART_OK else 'no'})")
                return _configure_solver(s, name)
        except Exception:
            continue
    raise RuntimeError(f"No se encontró solver en {SOLVER_ORDER}.")


import math as _math


def _solver_meta(m, res, wall_time, have_solution=False):
    """Extrae métricas portables entre backends de SOLVER_ORDER:
      · wall_time_s: tiempo de resolución (preferir el que reporte el solver;
        fallback universal = time.time() medido por _safe_solve).
      · bound: cota dual del MIP (upper_bound si max, lower_bound si min).
      · gap / gap_pct: brecha relativa incumbente-cota, calculada a mano
        (más fiable que el MIPGap específico de cada backend):
            gap = |incumbente − cota| / (|incumbente| + 1e-9).
    Si el solver no reporta cota (p.ej. GLPK legacy) ⇒ gap_pct='NA' explícito
    (distinto de 0 = "gap cerrado" y de vacío = "no medido")."""
    meta = dict(wall_time_s=float(wall_time), bound=None, gap=None, gap_pct="NA")
    if res is None:
        return meta
    # Tiempo reportado por el solver (si existe y es positivo) tiene prioridad.
    for attr in ("time", "wallclock_time", "wall_time", "user_time"):
        try:
            v = getattr(res.solver, attr, None)
            if v is not None and float(v) > 0.0:
                meta["wall_time_s"] = float(v)
                break
        except Exception:
            continue
    if not have_solution:
        return meta                              # sin incumbente cargado: sin gap
    # Incumbente = valor del objetivo activo con la solución ya cargada en m.
    try:
        inc = float(pyo.value(m.OBJ))
    except Exception:
        inc = None
    lb = ub = None
    try:
        lb = float(res.problem.lower_bound)
    except Exception:
        pass
    try:
        ub = float(res.problem.upper_bound)
    except Exception:
        pass

    def _fin(x):
        return x is not None and _math.isfinite(x)

    try:
        sense = m.OBJ.sense
    except Exception:
        sense = pyo.maximize
    bound = ub if sense == pyo.maximize else lb   # cota dual según el sentido
    if not _fin(bound):
        bound = lb if _fin(lb) else (ub if _fin(ub) else None)
    if _fin(bound):
        meta["bound"] = bound
        if inc is not None:
            meta["gap"] = abs(inc - bound) / (abs(inc) + 1e-9)
            meta["gap_pct"] = 100.0 * meta["gap"]
    return meta


def _safe_solve(solver, m, tee=False, warmstart=False):
    """Devuelve (status, res, meta) — meta con wall_time_s / bound / gap / gap_pct."""
    from pyomo.opt import TerminationCondition as _TC
    kw = dict(tee=tee, load_solutions=False)
    if warmstart and SOLVER_WARMSTART and _SOLVER_NAME in _WARMSTART_OK:
        kw["warmstart"] = True
    t0 = time.time()
    try:
        res = solver.solve(m, **kw)
    except Exception:
        try:                                     # reintento sin warmstart
            res = solver.solve(m, tee=tee, load_solutions=False)
        except Exception:
            return "error", None, _solver_meta(m, None, time.time() - t0)
    wall = time.time() - t0
    tc = res.solver.termination_condition
    if tc in (_TC.infeasible, _TC.infeasibleOrUnbounded, _TC.unbounded):
        return "infeasible", res, _solver_meta(m, res, wall)
    try:
        if len(res.solution) > 0:
            m.solutions.load_from(res)
            meta = _solver_meta(m, res, wall, have_solution=True)
            if tc == _TC.optimal:
                return "optimal", res, meta
            elif tc == _TC.maxTimeLimit:
                return "timelimit", res, meta
            else:
                return "feasible", res, meta
    except Exception:
        pass
    return "nosolution", res, _solver_meta(m, res, wall)


def set_objective(m, expr, sense):
    if m.component("OBJ") is not None:
        m.del_component("OBJ")
    m.OBJ = pyo.Objective(expr=expr, sense=sense)


def _open_sum(m):
    T_last = m._D["TT"][-1]
    return sum(m.y[k, T_last] for k in m._D["K"])


# ===========================================================================
# REPORTE Y PAYOFF
# ===========================================================================
def report(m, tag):
    D = m._D
    T_last = D["TT"][-1]
    if not D["K"]:
        return dict(tag=tag, n_open=0, opened=[], coverage_pct=0.0, NPV=0, Z1=0, Z2=0, Z3=0, phi=0, pi_individual={})

    opened = sorted(k for k in D["K"] if pyo.value(m.y[k, T_last]) > 0.5)
    coll = sum(pyo.value(m.b[i, j, mm, t])
               for (i, j) in [(i, j) for i, js in D["adj_ij"].items() for j in js]
               for mm in D["M"] for t in D["TT"])
    gen = sum(D["w"].get((i, mm, t), 0.0) for i in D["I"] for mm in D["M"] for t in D["TT"])
    pi_vals = {k: float(pyo.value(m.pi[k])) for k in opened}
    return dict(tag=tag, n_open=len(opened), opened=opened, coverage_pct=100.0 * coll / gen if gen else 0.0,
                NPV=pyo.value(m.NPV), Z1=pyo.value(m.Z1), Z2=pyo.value(m.Z2), Z3=pyo.value(m.Z3),
                phi=pyo.value(m.phi), pi_individual=pi_vals)


OBJ_SENSE = {"Z1": pyo.maximize, "Z2": pyo.minimize, "Z3": pyo.maximize,
             "Z4": pyo.maximize, "Z5": pyo.maximize}
OBJ_LABEL = {"Z1": "VAN peor planta", "Z2": "GEI neto", "Z3": "cobertura social",
             "Z4": "VAN total", "Z5": "output recuperado"}


def _obj_expr(m, key):
    return {"Z1": m.Z1, "Z2": m.Z2, "Z3": m.Z3, "Z4": m.NPV, "Z5": m.Z5}[key]


# Orden lexicográfico secundario FIJO, común a todas las filas de la tabla de
# pagos. Cada fila 'primary' se resuelve en el orden [primary] + (este orden sin
# primary), de modo que ECON_PRIMARY (Z4) encabeza el desempate salvo en su propia
# fila, donde ya ocupa su lugar natural (el primero). Usar el MISMO orden en todas
# las filas es lo que hace la tabla reproducible: las celdas fuera de la diagonal
# dependen del orden de desempate, y un orden distinto por fila daría rangos
# distintos para la misma instancia.
LEX_SECONDARY_ORDER = ["Z4", "Z2", "Z5", "Z1", "Z3"]


def _lex_row_order(primary, keys):
    """Orden de optimización de la fila 'primary' de la tabla de pagos."""
    return [primary] + [k for k in LEX_SECONDARY_ORDER if k != primary and k in keys]


def _add_lex_fix(m, name, key, z_prev):
    """Restricción temporal que mantiene el objetivo 'key' cerca de su óptimo
    z_prev mientras se optimiza el siguiente objetivo del orden lexicográfico.

    ⚠ Diferencia deliberada con `_set_eps`: allí la tolerancia se escala con R[k]
    (la amplitud del rango del objetivo k en la tabla de pagos), porque en
    `epsilon_constraint` la tabla YA está calculada. Aquí estamos justamente
    calculándola, así que R[k] no existe todavía: la tolerancia se escala con el
    propio óptimo, tol = SOLVER_MIP_GAP * max(1,|z_prev|). Es la misma idea —atar
    la holgura al gap real del solver— con la única escala disponible en esta
    etapa. Se usa DESIGUALDAD (no igualdad estricta) para no chocar con
    SOLVER_MIP_GAP: z_prev sólo está certificado dentro de esa banda, y una
    igualdad exacta puede volver infactible la fila entera."""
    tol = SOLVER_MIP_GAP * max(1.0, abs(z_prev))
    f = _obj_expr(m, key)
    con = (f >= z_prev - tol) if OBJ_SENSE[key] == pyo.maximize else (f <= z_prev + tol)
    m.add_component(name, pyo.Constraint(expr=con))


def lexicographic_payoff_table(m, solver, objs=("Z1", "Z2", "Z3", "Z4", "Z5")):
    """Tabla de pagos LEXICOGRÁFICA — Mavrotas (2009), Sec. 3.1.

    [M09] es explícito: la tabla de pagos debe construirse con optimización
    lexicográfica, no con óptimos individuales aislados. Optimizar cada objetivo
    por separado deja las celdas fuera de la diagonal en un punto arbitrario entre
    los muchos óptimos alternativos del MILP, y de ahí sale un nadir SOBRE-estimado
    (rangos inflados) o directamente un rango que no corresponde a ninguna solución
    eficiente. La grilla ε construida sobre esos rangos coloca puntos en zonas
    infactibles o dominadas.

    Cada fila 'primary' recorre `_lex_row_order(primary)`:
      1. Optimiza el primer objetivo del orden (= primary). z_prev := su óptimo.
      2. Para cada objetivo siguiente: añade una restricción temporal que amarra
         el objetivo ANTERIOR a z_prev (dentro del gap del solver), optimiza el
         siguiente, y actualiza z_prev.
      3. Registra en cells[primary] el valor de TODOS los objetivos en el punto
         final de la cadena (no sólo el primario).
      4. Borra las restricciones de la fila antes de pasar a la siguiente.

    El punto final de cada fila es un extremo (lexicográficamente) eficiente del
    frente: se registra como `lex_anchor_{primary}` y se siembra en `seed_pts` para
    garantizar que los vértices del frente aparezcan aunque la grilla ε los pierda
    por timeout.

    Coste: |objs| filas × hasta |objs| solves ≈ 20-25 solves.
    Devuelve (cells, ranges, anchors) — misma firma que la versión mono-objetivo.
    """
    keys = tuple(objs)
    cells = {}
    anchors = []
    t_all = time.time()
    first = True
    print(f"  payoff lexicográfico [M09 Sec. 3.1] — orden secundario: {LEX_SECONDARY_ORDER}")
    print("    " + " " * 9 + "".join(f"{k+'('+OBJ_LABEL.get(k,'')+')':>20s}" for k in keys))
    for primary in keys:
        order = _lex_row_order(primary, keys)
        fixes = []            # nombres de las restricciones temporales de esta fila
        z_prev = None
        meta = {"gap_pct": "NA", "bound": None, "wall_time_s": 0.0}
        for step, key in enumerate(order):
            if step > 0:
                # Amarra el objetivo anterior del orden a su óptimo recién hallado.
                name = f"lex_fix_{primary}_{step}"
                _add_lex_fix(m, name, order[step - 1], z_prev)
                fixes.append(name)
            set_objective(m, _obj_expr(m, key), OBJ_SENSE[key])
            cat, _, meta = _safe_solve(solver, m, tee=SOLVER_TEE, warmstart=not first)
            first = False
            if cat in ("infeasible", "nosolution", "error"):
                # El modelo conserva la solución del paso anterior (load_solutions=False
                # en el solve fallido): la fila se cierra ahí, con lo que ya se tiene.
                print(f"    fila {primary}: paso {step} ({key}) -> {cat}; "
                      f"se corta la cadena lexicográfica y se registra el último punto factible")
                break
            z_prev = float(pyo.value(_obj_expr(m, key)))

        achieved = {}
        for key in keys:
            try:
                achieved[key] = float(pyo.value(_obj_expr(m, key)))
            except Exception:
                achieved[key] = 0.0

        try:
            ap = _record_point(m, status=f"lex_anchor_{primary}",
                               gap_pct=meta["gap_pct"], bound=meta["bound"],
                               wall_time_s=meta["wall_time_s"])
            if ap["n_open"] > 0:
                anchors.append(ap)
        except Exception:
            pass

        for name in fixes:                        # limpieza: la fila siguiente parte limpia
            m.del_component(name)

        cells[primary] = {k: achieved[k] for k in keys}
        c = cells[primary]
        print(f"    lex {primary:5s}" + "".join(f"{c[k]:>20,.1f}" for k in keys))

    ranges = {k: (min(cells[o][k] for o in keys), max(cells[o][k] for o in keys)) for k in keys}
    print(f"\n  tabla de pagos completa en {time.time()-t_all:.1f}s")
    print("  rangos [min, max] por objetivo (entrada de la grilla ε):")
    for k in keys:
        lo, hi = ranges[k]
        amp = hi - lo
        # Mismo umbral que usa epsilon_constraint para partir active / loose_only.
        act = "activo" if amp > 1e-6 * max(1.0, abs(hi)) else "rango~0 (se fija flojo)"
        print(f"    {k} ({OBJ_LABEL.get(k,''):>16s}): [{lo:>16,.1f}, {hi:>16,.1f}]  "
              f"amplitud={amp:>16,.1f}  -> {act}")
    return cells, ranges, anchors


# ===========================================================================
# ε-CONSTRAINT (AUGMECON)
# ===========================================================================
def _ensure_augmecon(m, constrained):
    """Crea, por cada objetivo restringido, su RHS mutable eps_k, su holgura s_k≥0
    y la igualdad AUGMECON:  f_k − s_k = eps_k  (max, cota inferior f_k≥eps_k)  ó
    f_k + s_k = eps_k  (min, cota superior f_k≤eps_k). Idempotente.

    Réplica exacta de la ec. (3) de Mavrotas (2009): las ε-constraints se escriben
    como IGUALDADES con variable de holgura no negativa, no como desigualdades."""
    for key in constrained:
        if m.component(f"aug_s_{key}") is not None:
            continue
        big = 1e18 if OBJ_SENSE[key] == pyo.minimize else -1e18
        m.add_component(f"aug_eps_{key}", pyo.Param(initialize=big, mutable=True))
        m.add_component(f"aug_s_{key}", pyo.Var(domain=pyo.NonNegativeReals))
        f = _obj_expr(m, key)
        s = m.component(f"aug_s_{key}")
        eps = m.component(f"aug_eps_{key}")
        con = (f - s == eps) if OBJ_SENSE[key] == pyo.maximize else (f + s == eps)
        m.add_component(f"aug_c_{key}", pyo.Constraint(expr=con))


def _augmecon_ranges(rng, constrained):
    """R_k = amplitud (>=1) del rango del objetivo k en la tabla de pagos. Es el
    r_k de la ec. (6) de Mavrotas (2009), que normaliza cada holgura s_k en el
    término aumentado. También da la escala a la tolerancia ε de _set_eps."""
    return {k: max(rng[k][1] - rng[k][0], 1.0) for k in constrained}


def _set_eps(m, R, k, e):
    """Fija el RHS mutable eps_k de la ec. (3), aflojado una tolerancia numérica.

    ⚠ La tolerancia `ftol` NO forma parte de AUGMECON: la ec. (3) de Mavrotas
    (2009) fija el RHS exacto. Es un añadido de implementación exigido por
    resolver un MILP (y no un LP): las celdas de la tabla de pagos sólo están
    certificadas dentro de SOLVER_MIP_GAP, de modo que fijar eps_k exactamente en
    el óptimo individual de k puede volver infactible el extremo ideal. Se ata al
    gap real del solver — ver nota (b) en la cabecera de configuración."""
    ftol = EPS_FEAS_TOL_FRAC * SOLVER_MIP_GAP * R[k]
    v = (e - ftol) if OBJ_SENSE[k] == pyo.maximize else (e + ftol)
    m.component(f"aug_eps_{k}").set_value(float(v))


def _setup_augmecon_objective(m, primary, constrained, R, augw=ECON_AUGW):
    """Arma el objetivo AUGMECON con el sentido del primario. Compartido por
    epsilon_constraint y --time-curve para garantizar que la curva de optimalidad
    resuelve EL MISMO subproblema.

    Ec. (6) de Mavrotas (2009):   max ( f1(x) + eps * Σ_k s_k / r_k )
    con r_k = amplitud del rango de f_k en la tabla de pagos (aquí R[k]) y
    eps = augw ∈ [1e-3, 1e-6]. Para un primario a minimizar se usa la forma
    equivalente  min ( f1(x) − eps * Σ_k s_k / r_k ).

    TIEBREAK_OPEN_W (por defecto 0.0) añadiría un tercer término que NO existe en
    la ec. (6). Se mantiene apagado para que el objetivo sea textualmente el del
    paper; ver nota (c) en la cabecera de configuración."""
    prim = _obj_expr(m, primary)
    psense = OBJ_SENSE[primary]
    aug = augw * sum(m.component(f"aug_s_{k}") / R[k] for k in constrained)
    expr = (prim + aug) if psense == pyo.maximize else (prim - aug)
    if TIEBREAK_OPEN_W:
        # Fuera de AUGMECON: penalización explícita de desempate por nº de plantas.
        tb = TIEBREAK_OPEN_W * _open_sum(m)
        expr = (expr - tb) if psense == pyo.maximize else (expr + tb)
    set_objective(m, expr, psense)
    return psense


def _set_time_limit(solver, tl):
    """Reajusta el límite de tiempo del solver activo entre solves (para la curva
    de optimalidad, que barre TIME_LIMITS crecientes sobre la misma celda)."""
    name = _SOLVER_NAME
    tl = float(tl)
    if name == "gurobi":
        solver.options["TimeLimit"] = tl
    elif name == "cbc":
        solver.options["seconds"] = tl
    elif name == "appsi_highs":
        solver.options["time_limit"] = tl
    elif name == "glpk":
        solver.options["tmlim"] = int(tl)


def _record_point(m, **extra):
    D = m._D
    I, J, K, MM, TT = D["I"], D["J"], D["K"], D["M"], D["TT"]
    T_last = TT[-1]
    price, cprod, eprod = D["price"], D["cprod"], D["eprod"]
    ctr, eptr = D["ctr"], D["eptr"]
    djk, akj = D["A_JK"], D["adj_kj"]

    yv = m.y.extract_values()
    xv = m.x.extract_values()
    qv = m.q.extract_values()
    bv = m.b.extract_values()

    def _is_open(k, t):
        return 1 if (yv.get((k, t)) or 0.0) > 0.5 else 0

    opened = sorted(k for k in K if _is_open(k, T_last))
    open_period = {}
    for k in opened:
        for t in TT:
            if _is_open(k, t):
                open_period[k] = int(t)
                break

    temporal_y = {str(k): {str(t): _is_open(k, t) for t in TT} for k in opened}

    temporal_econ = {}
    for k in opened:
        rows = {}
        for t in TT:
            q_m = {mm: float(qv.get((k, mm, t)) or 0.0) for mm in MM}
            x_in = {mm: sum(float(xv.get((j, k, mm, t)) or 0.0) for j in akj[k]) for mm in MM}
            rev = sum(price[mm] * q_m[mm] for mm in MM)
            tr = ctr * 1e-3 * sum(djk[(j, k)] * float(xv.get((j, k, mm, t)) or 0.0) for j in akj[k] for mm in MM)
            pr = sum(cprod[mm] * x_in[mm] for mm in MM)
            emis = (eptr * 1e-3 * sum(djk[(j, k)] * float(xv.get((j, k, mm, t)) or 0.0) for j in akj[k] for mm in MM)
                    + sum(eprod[mm] * x_in[mm] for mm in MM))
            if LCA_CREDIT:
                emis -= sum(LCA_CREDIT_M.get(mm, 0.0) * q_m[mm] for mm in MM)
            if USE_DISPOSAL:                       # crédito por vertedero evitado al recuperar
                emis -= sum(DISPOSAL_EF_M.get(mm, 0.0) * q_m[mm] for mm in MM)
            rows[str(t)] = dict(
                y=_is_open(k, t),
                q_PET=q_m.get(1, 0.0), q_HDPE=q_m.get(2, 0.0), q_PP=q_m.get(3, 0.0),
                production=sum(q_m.values()),
                revenue=rev, transport=tr, processing=pr,
                margin=rev - tr - pr, emissions=emis)
        temporal_econ[str(k)] = rows

    pi_individual = {k: float(pyo.value(m.pi[k])) for k in opened}
    worst_pi = min(pi_individual.values()) if pi_individual else 0.0

    hub_timeline = {}
    zone_collection = {} if EXPORT_ZONE_DETAIL else None
    for (i, j, mm, t), v in bv.items():
        if not v or v <= 1e-9:
            continue
        fv = float(v)
        hj = hub_timeline.setdefault(str(j), {})
        hj[str(t)] = hj.get(str(t), 0.0) + fv
        if zone_collection is not None:
            zone_collection[i] = zone_collection.get(i, 0.0) + fv

    fjk = {}
    for (j, k, mm, t), v in xv.items():
        if v and v > 1e-9:
            fjk[(j, k)] = fjk.get((j, k), 0.0) + float(v)
    flows_jk = [[int(j), int(k), float(f)] for (j, k), f in sorted(fjk.items())]

    collected_total = sum(v for tl in hub_timeline.values() for v in tl.values())
    gen_total = sum(D["Gamma"].values())

    pt = dict(
        Z1=float(pyo.value(m.Z1)), Z2=float(pyo.value(m.Z2)), Z3=float(pyo.value(m.Z3)), Z5=float(pyo.value(m.Z5)),
        phi=float(pyo.value(m.phi)), pi_individual=pi_individual, worst_pi=worst_pi,
        coverage=float(pyo.value(m.coverage)), NPV=float(pyo.value(m.NPV)),
        Z4=float(pyo.value(m.NPV)),            # VAN total (4º objetivo) = alias de NPV
        collected_total=collected_total,
        coverage_pct=(100.0 * collected_total / gen_total if gen_total else 0.0),
        n_open=len(opened), opened=opened,
        open_period=open_period, temporal_y=temporal_y, temporal_econ=temporal_econ,
        hub_timeline=hub_timeline, flows_jk=flows_jk,
        _zone_collection=zone_collection,
    )
    pt.update(extra)
    return pt


def _pareto_filter(pts, tol=1e-6, objs=None):
    """Filtro de no-dominancia post-hoc sobre los puntos obtenidos.

    ⚠ Salvaguarda exigida por trabajar con MILP + gap, NO parte de AUGMECON. La
    Proposición de Mavrotas (2009), Sec. 3.2, garantiza que toda solución óptima
    de (3) es eficiente, pero su demostración supone OPTIMALIDAD EXACTA. Con
    SOLVER_MIP_GAP = 0.005 cada celda se resuelve sólo dentro de un 0.5% del
    óptimo, así que la garantía no aplica de forma estricta: dos celdas vecinas
    pueden devolver incumbentes mutuamente dominados dentro de esa banda. Este
    filtro descarta las dominadas a posteriori y acota la desviación. Debe
    declararse como tal en la metodología del manuscrito."""
    if not pts:
        return []
    objs = [k for k in (objs or ([ECON_PRIMARY] + list(ECON_CONSTRAINED_ALL))) if k in pts[0]]
    uniq, seen = [], set()
    for p in pts:
        key = tuple(round(p[k], 3) for k in objs)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(p)

    def _be(a, b, k):   # a al menos tan bueno como b en k (según su sentido)
        return (a[k] >= b[k] - tol) if OBJ_SENSE[k] == pyo.maximize else (a[k] <= b[k] + tol)

    def _st(a, b, k):   # a estrictamente mejor que b en k
        return (a[k] > b[k] + tol) if OBJ_SENSE[k] == pyo.maximize else (a[k] < b[k] - tol)

    def dom(a, b):
        return all(_be(a, b, k) for k in objs) and any(_st(a, b, k) for k in objs)
    return [p for p in uniq if not any(dom(q, p) for q in uniq if q is not p)]


def epsilon_constraint(m, solver, cells, rng, primary=None, constrained=None,
                       npts=ECON_NPTS, npts_map=None, augw=ECON_AUGW, seed_pts=None):
    """AUGMECON de p objetivos: 'primary' se optimiza (ec. (6)); el resto entran como
    ε-constraints aumentadas (ec. (3)). Bucles anidados (producto de los objetivos
    externos) × bucle interno = último de 'constrained'.

    `npts_map` da el q_i de cada objetivo ([M09] Sec. 4): la resolución de la grilla
    es POR OBJETIVO, no uniforme. Los objetivos ausentes del mapa usan `npts`.

    Fidelidad al método (ver cabecera de configuración):
      · Mavrotas (2009), Sec. 3.3 / Fig. 4 — early-exit al volverse INFACTIBLE el
        subproblema: el resto de la fila del bucle interno también lo será.
      · Mavrotas & Florios (2013) — bypass de puntos de grilla FACTIBLES según la
        magnitud de la holgura interna. Es AUGMECON2, no AUGMECON: sólo se aplica
        si AUGMECON2_BYPASS=True, y en ese caso el manuscrito debe citarlo."""
    import itertools
    primary = primary or ECON_PRIMARY
    constrained = list(constrained or ECON_CONSTRAINED_ALL)
    _ensure_augmecon(m, constrained)

    # Rangos del payoff. Los objetivos con rango informativo se grillan (active); los
    # de rango ~nulo (sin trade-off) se fijan en su extremo flojo/nadir (no vinculante)
    # y no consumen dimensiones de grilla — evita corners infactibles en instancias
    # con frente casi puntual.
    R = _augmecon_ranges(rng, constrained)
    active = [k for k in constrained if (rng[k][1] - rng[k][0]) > 1e-6 * max(1.0, abs(rng[k][1]))]
    loose_only = [k for k in constrained if k not in active]

    _setup_augmecon_objective(m, primary, constrained, R, augw=augw)

    def _rel(k):
        return "≤" if OBJ_SENSE[k] == pyo.minimize else "≥"

    for k in loose_only:                          # fija degenerados en su nadir (flojo)
        lo, hi = rng[k]
        _set_eps(m, R, k, lo if OBJ_SENSE[k] == pyo.maximize else hi)

    # q_i por objetivo ([M09] Sec. 4). Los ausentes del mapa caen en 'npts'.
    npts_map = dict(npts_map or {})
    q = {k: max(1, int(npts_map.get(k, npts))) for k in active}

    # Dirección de la grilla, Mavrotas (2009) Sec. 3.3: para un objetivo a maximizar
    # se parte del mínimo del rango y se sube (la ε-constraint se va apretando);
    # para uno a minimizar, al revés.
    step, grid = {}, {}
    for k in active:
        lo, hi = rng[k]
        step[k] = (hi - lo) / (q[k] - 1) if q[k] > 1 else 0.0
        grid[k] = ([lo + i * step[k] for i in range(q[k])] if OBJ_SENSE[k] == pyo.maximize
                   else [hi - i * step[k] for i in range(q[k])])

    n_total = 1
    for k in active:
        n_total *= q[k]
    method = ("AUGMECON2 [Mavrotas & Florios 2013]" if AUGMECON2_BYPASS
              else "AUGMECON [Mavrotas 2009]")
    accel = "early-exit + bypass + warmstart" if AUGMECON2_BYPASS else "early-exit + warmstart"
    qtxt = " × ".join(f"{k}:{q[k]}" for k in active) if active else "—"
    print(f"\n== {method}: primary={primary}, ε-activos={active}"
          + (f", fijos(rango~0)={loose_only}" if loose_only else "")
          + f", {n_total} solves máx ({accel}) ==")
    print(f"  densidad de grilla q_i [M09 Sec. 4]: {qtxt}  ->  n_total = {n_total}")

    raw, k_done, first = [], 0, True
    if seed_pts:                                  # anclas del payoff (extremos garantizados)
        raw.extend(seed_pts)
        print(f"  sembradas {len(seed_pts)} anclas del payoff en el frente")

    if not active:
        # Frente de un solo punto: un único solve con todos los objetivos flojos.
        cat, _, meta = _safe_solve(solver, m, tee=SOLVER_TEE, warmstart=False)
        if cat not in ("infeasible", "nosolution", "error"):
            pt = _record_point(m, status=cat, gap_pct=meta["gap_pct"],
                               bound=meta["bound"], wall_time_s=meta["wall_time_s"])
            print(f"  [1/1] punto único  Z1={pt['Z1']:,.0f} Z2={pt['Z2']:,.0f} "
                  f"Z3={pt['Z3']:.1f} Z4={pt.get('Z4', pt['NPV']):,.0f} open={pt['n_open']}")
            raw.append(pt)
        else:
            print(f"  [1/1] -> {cat}")
    else:
        inner = active[-1]
        outer = active[:-1]
        istep = step[inner]
        q_in = q[inner]                           # q_i del bucle interno
        outer_combos = list(itertools.product(*[grid[k] for k in outer])) if outer else [()]
        for combo in outer_combos:
            for k, e in zip(outer, combo):
                _set_eps(m, R, k, e)
            oc = "  ".join(f"{k}{_rel(k)}{e:>10,.0f}" for k, e in zip(outer, combo))
            j = 0
            while j < q_in:
                e_in = grid[inner][j]
                k_done += 1
                _set_eps(m, R, inner, e_in)
                cat, _, meta = _safe_solve(solver, m, tee=SOLVER_TEE, warmstart=not first)
                first = False
                tag = f"  [{k_done:>4d}/{n_total}] {oc}  {inner}{_rel(inner)}{e_in:>10,.1f}"
                if cat == "infeasible":
                    # Early-exit de Mavrotas (2009), Sec. 3.3 y Fig. 4: al volverse
                    # infactible, los ε internos restantes (más apretados) también lo son.
                    skip = q_in - j - 1
                    print(f"{tag}  -> infeasible (early-exit, salta {skip})")
                    k_done += skip
                    break
                if cat in ("nosolution", "error"):
                    print(f"{tag}  -> {cat} (sin incumbente en {SOLVER_TIME_LIMIT}s; no descarta la fila)")
                    j += 1
                    continue
                extra = {f"eps_{k}": float(e) for k, e in zip(outer, combo)}
                extra[f"eps_{inner}"] = float(e_in)
                extra["status"] = cat
                extra["gap_pct"] = meta["gap_pct"]
                extra["bound"] = meta["bound"]
                extra["wall_time_s"] = meta["wall_time_s"]
                pt = _record_point(m, **extra)
                # Bypass de Mavrotas & Florios (2013) — AUGMECON2, NO Mavrotas (2009):
                # la holgura del objetivo interno cubre bp = floor(s_in/step_in) puntos
                # de grilla siguientes con la MISMA solución ⇒ se omiten. Desactivable
                # con AUGMECON2_BYPASS=False para caer en AUGMECON puro [M09].
                bp = 0
                if AUGMECON2_BYPASS and istep > 1e-12:
                    try:
                        s_in = float(pyo.value(m.component(f"aug_s_{inner}")))
                    except Exception:
                        s_in = 0.0
                    bp = max(0, min(int(s_in / istep), q_in - 1 - j))
                note = f"  [bypass +{bp}]" if bp else ""
                gp = meta["gap_pct"]
                gtxt = f"{gp:>5.2f}%" if isinstance(gp, (int, float)) else f"{gp:>6}"
                flag = "  <TL>" if cat == "timelimit" else ""
                print(f"{tag}  Z1={pt['Z1']:>11,.0f} Z2={pt['Z2']:>11,.0f} Z3={pt['Z3']:>8.1f} "
                      f"Z4={pt.get('Z4', pt['NPV']):>11,.0f} open={pt['n_open']} "
                      f"gap={gtxt} t={meta['wall_time_s']:>5.0f}s{flag}{note}")
                raw.append(pt)
                k_done += bp
                j += bp + 1

    # La no-dominancia se evalúa en el espacio de objetivos que REALMENTE varía
    # (primario + activos). Incluir los de rango~0 sólo debilitaría el filtro: son
    # casi constantes, así que nunca desempatan y su ruido numérico haría pasar
    # puntos dominados. Se siguen calculando y reportando en cada punto.
    front = _pareto_filter(raw, objs=[primary] + active)
    print(f"\n  {len(raw)} solves con solución -> {len(front)} puntos no dominados "
          f"(dominancia sobre {[primary] + active})")
    return front, raw, cells


# ===========================================================================
# AUDITORÍA DE GAP / CURVA DE OPTIMALIDAD
# ===========================================================================
def _is_problematic(p, tol_pct=GAP_AUDIT_TOL_PCT):
    """Una celda es problemática solo si su gap medido supera tol_pct."""
    gp = p.get("gap_pct")
    return (
        isinstance(gp, (int, float))
        and not isinstance(gp, bool)
        and gp > tol_pct
    )


def _present_eps_keys(raw, constrained):
    """Objetivos restringidos que REALMENTE tienen ε registrado en los puntos =
    los 'activos' grillados (los 'loose' de rango~0 se fijan y no se registran)."""
    for p in raw:
        present = [k for k in constrained if f"eps_{k}" in p]
        if present:
            return present
    return []


def select_problematic_cells(raw, constrained=None, tol_pct=GAP_AUDIT_TOL_PCT):
    """Extrae de 'raw' las celdas de grilla (con eps_{k} conocidos) que no cerraron
    gap, deduplicadas por su tupla de ε. Devuelve la lista tal cual (cada elemento
    conserva eps_{k}, status, gap_pct) — entrada exacta del modo --time-curve."""
    constrained = list(constrained or ECON_CONSTRAINED_ALL)
    eps_keys = _present_eps_keys(raw, constrained)
    cells, seen = [], set()
    for p in raw:
        if not eps_keys or not all(f"eps_{k}" in p for k in eps_keys):
            continue                              # anclas del payoff: sin ε, se omiten
        if not _is_problematic(p, tol_pct):
            continue
        key = tuple(round(float(p[f"eps_{k}"]), 3) for k in eps_keys)
        if key in seen:
            continue
        seen.add(key)
        cells.append(p)
    return cells


def summarize_gaps(raw, constrained=None, tol_pct=GAP_AUDIT_TOL_PCT):
    """Imprime el resumen de gap por celda de grilla y lista las problemáticas
    (entrada del Fix 4). Devuelve la lista de celdas problemáticas."""
    constrained = list(constrained or ECON_CONSTRAINED_ALL)
    eps_keys = _present_eps_keys(raw, constrained)
    grid_pts = [p for p in raw if eps_keys and all(f"eps_{k}" in p for k in eps_keys)]
    print("\n== Auditoría de gap por celda de grilla ==")
    if not grid_pts:
        print("  (sin celdas de grilla con ε registrados)")
        return []
    n_tl = sum(1 for p in grid_pts if p.get("status") == "timelimit")
    n_na = sum(1 for p in grid_pts if not isinstance(p.get("gap_pct"), (int, float)))
    gaps = [p["gap_pct"] for p in grid_pts if isinstance(p.get("gap_pct"), (int, float))]
    if gaps:
        print(f"  celdas={len(grid_pts)}  gap medido: min={min(gaps):.3f}%  "
              f"máx={max(gaps):.3f}%  medio={sum(gaps)/len(gaps):.3f}%  "
              f"(timelimit={n_tl}, gap='NA'={n_na})")
    else:
        print(f"  celdas={len(grid_pts)}  sin gap medido (timelimit={n_tl}, NA={n_na})")
    prob = select_problematic_cells(raw, constrained, tol_pct)
    print(f"\n  celdas problemáticas (timelimit o gap>{tol_pct}%): {len(prob)}")
    for i, p in enumerate(prob, 1):
        eps_txt = "  ".join(f"eps_{k}={float(p[f'eps_{k}']):,.1f}" for k in eps_keys)
        gp = p.get("gap_pct")
        gtxt = f"{gp:.3f}%" if isinstance(gp, (int, float)) else str(gp)
        print(f"    [{i}] {eps_txt}  status={p.get('status')}  gap={gtxt}  "
              f"t={p.get('wall_time_s', '?')}s  Z4={p.get('Z4', p.get('NPV', 0)):,.0f}")
    if not prob:
        print("    (todas las celdas cerraron gap: no hace falta --time-curve)")
    return prob


def run_time_curve(m, solver, rng, cells, outdir, primary=None, constrained=None,
                   time_limits=None, stop_pct=TIME_CURVE_GAP_STOP_PCT):
    """Curva tiempo-vs-gap. Para cada celda problemática resuelve EL MISMO
    subproblema AUGMECON (mismos ε, misma función objetivo) con límites de tiempo
    crecientes, reusando el incumbente previo como warmstart dentro de la celda.
    Corte temprano si el gap baja de stop_pct. Escribe optimality_curve.csv."""
    primary = primary or ECON_PRIMARY
    constrained = list(constrained or ECON_CONSTRAINED_ALL)
    time_limits = list(time_limits or TIME_CURVE_LIMITS)
    os.makedirs(outdir, exist_ok=True)

    _ensure_augmecon(m, constrained)
    R = _augmecon_ranges(rng, constrained)
    active = [k for k in constrained if (rng[k][1] - rng[k][0]) > 1e-6 * max(1.0, abs(rng[k][1]))]
    loose = [k for k in constrained if k not in active]
    _setup_augmecon_objective(m, primary, constrained, R)
    for k in loose:                               # degenerados fijos en su nadir (igual que la grilla)
        lo, hi = rng[k]
        _set_eps(m, R, k, lo if OBJ_SENSE[k] == pyo.maximize else hi)

    print(f"\n== Modo curva de optimalidad: {len(cells)} celdas × {len(time_limits)} límites "
          f"(early-stop gap<{stop_pct}%) ==")
    print(f"   ε-activos={active}" + (f", fijos={loose}" if loose else ""))
    rows = []
    for ci, cell in enumerate(cells, 1):
        eps_cell = {k: float(cell[f"eps_{k}"]) for k in active if f"eps_{k}" in cell}
        for k in eps_cell:
            _set_eps(m, R, k, eps_cell[k])
        eps_hdr = "  ".join(f"{k}={eps_cell[k]:,.1f}" for k in eps_cell)
        print(f"\n  celda [{ci}/{len(cells)}]  {eps_hdr}")
        warm = False
        for tl in time_limits:
            _set_time_limit(solver, tl)
            cat, res, meta = _safe_solve(solver, m, tee=SOLVER_TEE, warmstart=warm)
            gp = meta["gap_pct"]
            gtxt = f"{gp:.3f}%" if isinstance(gp, (int, float)) else str(gp)
            if cat in ("infeasible", "nosolution", "error"):
                print(f"    TL={tl:>5}s -> {cat}  (t={meta['wall_time_s']:.0f}s)")
                rows.append([ci] + [eps_cell.get(k, "") for k in active]
                            + [tl, cat, gp, meta["wall_time_s"], "", "", "", "", ""])
                if cat == "infeasible":
                    break                         # ε fijos: infeasible no mejora con más tiempo
                continue
            rep = report(m, f"tc_{ci}_{tl}")
            rows.append([ci] + [eps_cell.get(k, "") for k in active]
                        + [tl, cat, gp, meta["wall_time_s"],
                           rep["NPV"], rep["Z1"], rep["Z2"], rep["Z3"], rep["n_open"]])
            print(f"    TL={tl:>5}s -> {cat:<9} gap={gtxt:<9} t={meta['wall_time_s']:>5.0f}s "
                  f"NPV={rep['NPV']:>12,.0f} Z2={rep['Z2']:>11,.0f} Z3={rep['Z3']:>8.1f} "
                  f"open={rep['n_open']}")
            warm = True                           # reusar incumbente en el siguiente TL de la celda
            if isinstance(gp, (int, float)) and gp < stop_pct:
                print(f"    -> gap cerrado (<{stop_pct}%) en {tl}s; salto los TL restantes")
                break

    hdr = (["cell"] + [f"eps_{k}" for k in active]
           + ["time_limit_s", "termination", "gap_pct", "wall_time_s",
              "NPV", "Z1", "Z2", "Z3", "n_open"])
    path = os.path.join(outdir, "optimality_curve.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(hdr)
        wr.writerows(rows)
    print(f"\n  curva de optimalidad → {path}  ({len(rows)} filas)")
    if MAKE_PLOTS:
        plot_optimality_curve(rows, active, os.path.join(outdir, "optimality_curve.png"))
    return rows


def plot_optimality_curve(rows, active, outfile):
    """gap_pct vs wall_time_s, una serie (con línea) por celda. A diferencia del
    frente de Pareto, aquí SÍ se conectan los puntos: es una curva de convergencia."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print("  matplotlib no disponible; omito optimality_curve.png")
        return
    n_eps = len(active)
    series = {}
    for r in rows:
        ci = r[0]
        gp = r[2 + n_eps + 2]                     # cell, eps..., time_limit, termination, gap_pct
        t = r[2 + n_eps + 3]                      # wall_time_s
        if not isinstance(gp, (int, float)) or not isinstance(t, (int, float)):
            continue
        series.setdefault(ci, []).append((t, gp))
    if not series:
        print("  (sin gaps numéricos que graficar en la curva de optimalidad)")
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    for ci, pts in sorted(series.items()):
        pts = sorted(pts)
        ax.plot([p[0] for p in pts], [p[1] for p in pts],
                marker="o", linewidth=1.4, markersize=5, label=f"celda {ci}")
    ax.set_xlabel("tiempo de resolución (s)")
    ax.set_ylabel("gap relativo (%)")
    ax.set_title("Curva de optimalidad: gap vs tiempo por celda")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140)
    plt.close(fig)
    print(f"  figura curva → {outfile}")


# ===========================================================================
# GUARDAR / CARGAR
# ===========================================================================
def _supply_seasonality(D):
    by_month = {}
    for (i, mm, t), w in D["w"].items():
        mo = (t - 1) % 12 + 1
        by_month[mo] = by_month.get(mo, 0.0) + w
    months = sorted(by_month)
    if not months:
        return [1.0] * 12
    vals = [by_month[mo] for mo in months]
    mu = sum(vals) / len(vals)
    return [v / mu for v in vals] if mu else [1.0] * len(vals)


def _build_geo(D):
    return dict(
        coords_i={str(i): list(D["coord_i"][i]) for i in D.get("coord_i", {})},
        coords_j={str(j): list(D["coord_j"][j]) for j in D.get("coord_j", {})},
        coords_k={str(k): list(D["coord_k"][k]) for k in D.get("coord_k", {})},
        A_JK=[[j, k] for (j, k) in D.get("A_JK", {})],
        seasonality=_supply_seasonality(D),
    )


def _strip_heavy(pt, keep_full):
    out = {k: v for k, v in pt.items() if k != "_zone_collection"}
    if not keep_full:
        for k in ("hub_timeline", "flows_jk"):
            out.pop(k, None)
    return out


def save_frontier(front, raw, cells, outdir, geo=None):
    os.makedirs(outdir, exist_ok=True)
    front_ids = {id(p) for p in front}
    payload = dict(
        params=dict(ZONE=ZONE, SCENARIO=SCENARIO, HORIZON=HORIZON, ALPHA=ALPHA, D_BAR_CENTER_K=D_BAR_CENTER_K,
                    D_BAR_IJ=D_BAR_IJ, D_BAR_JK=D_BAR_JK, D_BAR_KK_KM=D_BAR_KK_KM,
                    N_CAP=N_CAP, BETA_MIN=BETA_MIN, WASTE_FACTOR=WASTE_FACTOR,
                    LCA_CREDIT=LCA_CREDIT, USE_PRODUCT_PRICES=USE_PRODUCT_PRICES,
                    YIELD_MODE=YIELD_MODE, YIELD_SCALE=YIELD_SCALE, ETA_OVERRIDE=ETA_OVERRIDE,
                    USE_DISPOSAL=USE_DISPOSAL, DISPOSAL_EF_M=DISPOSAL_EF_M, EPS_TR_SCALE=EPS_TR_SCALE,
                    SOLVER_TIME_LIMIT=SOLVER_TIME_LIMIT, SOLVER_MIP_GAP=SOLVER_MIP_GAP,
                    ECON_PRIMARY=ECON_PRIMARY, ECON_CONSTRAINED_ALL=ECON_CONSTRAINED_ALL,
                    ECON_NPTS=ECON_NPTS, ECON_NPTS_MAP=ECON_NPTS_MAP,
                    LEX_SECONDARY_ORDER=LEX_SECONDARY_ORDER,
                    # Fidelidad al método (ver cabecera): qué se desvía de Mavrotas (2009).
                    ECON_AUGW=ECON_AUGW, AUGMECON2_BYPASS=AUGMECON2_BYPASS,
                    EPS_FEAS_TOL_FRAC=EPS_FEAS_TOL_FRAC, TIEBREAK_OPEN_W=TIEBREAK_OPEN_W),
        payoff=cells,
        pareto=[_strip_heavy(p, keep_full=True) for p in front],
        raw=[_strip_heavy(p, keep_full=(id(p) in front_ids)) for p in raw],
        geo=geo or {},
    )
    path = os.path.join(outdir, "econ_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    print(f"  resultados → {path}")
    return path


def save_detailed_csv(front, D, outdir):
    cdir = os.path.join(outdir, "detalle_csv")
    os.makedirs(cdir, exist_ok=True)
    coord_k = D.get("coord_k", {})
    coord_j = D.get("coord_j", {})
    coord_i = D.get("coord_i", {})
    Gamma = D.get("Gamma", {})

    def _w(name, header, rows):
        with open(os.path.join(cdir, name), "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(header)
            wr.writerows(rows)

    # Columnas ε dinámicas: se leen las claves reales eps_{k} que escribe
    # epsilon_constraint (antes se leía eps2/eps3, que NUNCA existían → vacías).
    # Iterar sobre ECON_CONSTRAINED_ALL evita que se rompa si el set de objetivos
    # restringidos cambia.
    eps_keys = list(ECON_CONSTRAINED_ALL)
    rows = []
    for s, p in enumerate(front, 1):
        rows.append([s, p["Z1"], p["Z2"], p["Z3"], p["NPV"], p.get("Z5", 0.0),
                     p.get("coverage_pct", 0.0), p.get("collected_total", 0.0),
                     p["n_open"], " ".join(map(str, p["opened"]))]
                    + [p.get(f"eps_{k}", "") for k in eps_keys]
                    + [p.get("gap_pct", "NA"), p.get("bound", ""),
                       p.get("wall_time_s", ""), p.get("status", "")])
    _w("soluciones.csv",
       ["sol", "Z1_VAN_EUR", "Z2_GHG_kgCO2e", "Z3_cobertura", "Z4_NPV_EUR", "Z5_output_kg",
        "cobertura_pct", "colectado_total_kg", "n_plantas", "plantas"]
       + [f"eps_{k}" for k in eps_keys]
       + ["gap_pct", "bound", "wall_time_s", "estado"], rows)

    rows = []
    for s, p in enumerate(front, 1):
        for k in p["opened"]:
            x, y, nm = coord_k.get(int(k), (None, None, ""))
            prod = sum(r["production"] for r in p["temporal_econ"].get(str(k), {}).values())
            rows.append([s, k, nm, x, y, p["open_period"].get(int(k), ""),
                         p["pi_individual"].get(int(k), p["pi_individual"].get(str(k), "")),
                         prod])
    _w("plantas_apertura.csv",
       ["sol", "k", "nombre", "x", "y", "apertura_t", "VAN_planta_EUR", "produccion_total_kg"], rows)

    rows = []
    for s, p in enumerate(front, 1):
        for k, serie in p["temporal_econ"].items():
            for t, r in serie.items():
                rows.append([s, int(k), int(t), r["y"], r["q_PET"], r["q_HDPE"], r["q_PP"],
                             r["production"], r["revenue"], r["transport"], r["processing"],
                             r["margin"], r["emissions"]])
    _w("plantas_temporal.csv",
       ["sol", "k", "t", "abierta", "q_PET_kg", "q_HDPE_kg", "q_PP_kg", "produccion_kg",
        "ingreso_EUR", "transporte_EUR", "procesamiento_EUR", "margen_EUR", "emisiones_kgCO2e"], rows)

    rows = []
    for s, p in enumerate(front, 1):
        for j, serie in p.get("hub_timeline", {}).items():
            x, y = coord_j.get(int(j), (None, None))
            for t, v in serie.items():
                rows.append([s, int(j), x, y, int(t), v])
    _w("hubs_temporal.csv", ["sol", "j", "x", "y", "t", "influjo_kg"], rows)

    rows = []
    for s, p in enumerate(front, 1):
        for j, k, f in p.get("flows_jk", []):
            xj, yj = coord_j.get(int(j), (None, None))
            xk, yk, nk = coord_k.get(int(k), (None, None, ""))
            rows.append([s, int(j), int(k), f, xj, yj, xk, yk])
    _w("flujos_jk.csv", ["sol", "j", "k", "flujo_kg", "x_j", "y_j", "x_k", "y_k"], rows)

    if EXPORT_ZONE_DETAIL and any(p.get("_zone_collection") for p in front):
        rows = []
        for s, p in enumerate(front, 1):
            zc = p.get("_zone_collection") or {}
            for i, gen in Gamma.items():
                x, y, _w_i = coord_i.get(int(i), (None, None, None))
                col = zc.get(i, 0.0)
                rows.append([s, int(i), x, y, gen, col, (col / gen if gen else 0.0)])
        _w("zonas_cobertura.csv",
           ["sol", "i", "x", "y", "generado_kg", "colectado_kg", "cobertura_frac"], rows)

    print(f"  CSV detallados → {cdir}/")


def load_frontier(path=None):
    if path is None:
        path = os.path.join(ECON_OUTDIR, "econ_results.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"No se encontró: {path}\nCorre primero sin --plot.")
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    front = data["pareto"]
    raw = data.get("raw", [])
    geo = data.get("geo", {})
    params = data.get("params", {})
    print(f"  cargado: {len(front)} puntos Pareto, {len(raw)} solves  ({path})")
    if params:
        print("  params: " + "  ".join(f"{k}={v}" for k, v in params.items()))
    return front, raw, geo, params


# ===========================================================================
# VISUALIZACIÓN
# ===========================================================================
_UNITS = {"Z1": "EUR", "Z2": "kg CO2eq", "Z3": "", "Z4": "EUR", "Z5": "kg"}


def _minmax(vals):
    a = np.asarray(vals, dtype=float)
    lo, hi = float(a.min()), float(a.max())
    normed = np.full_like(a, 0.5) if hi - lo < 1e-9 else (a - lo) / (hi - lo)
    return lo, hi, normed


def _scale(vals, lo, hi):
    a = np.asarray(vals, dtype=float)
    return np.full_like(a, 0.5) if hi - lo < 1e-9 else (a - lo) / (hi - lo)


def _axis_label(key, lo, hi):
    unit = _UNITS.get(key, "")
    fmt = ",.0f" if abs(hi) > 100 else ".4f"
    suffix = f" {unit}" if unit else ""
    return f"{key} (norm)  [0={lo:{fmt}}{suffix} | 1={hi:{fmt}}{suffix}]"


def plot_2d(front, xk, yk, outfile, title, raw=None):
    import matplotlib.pyplot as plt
    basis = raw if raw else front
    xlo, xhi, _ = _minmax([p[xk] for p in basis])
    ylo, yhi, _ = _minmax([p[yk] for p in basis])
    pts = sorted(front, key=lambda p: p[xk])

    fig, ax = plt.subplots(figsize=(7, 5))
    if raw:
        ax.scatter(_scale([p[xk] for p in raw], xlo, xhi),
                   _scale([p[yk] for p in raw], ylo, yhi),
                   s=18, color="lightgrey", edgecolors="none", zorder=1)
    fx = _scale([p[xk] for p in pts], xlo, xhi)
    fy = _scale([p[yk] for p in pts], ylo, yhi)
    ax.scatter(fx, fy, s=60, color="steelblue", edgecolors="k", linewidths=0.6, zorder=3)
    for i, (x, y) in enumerate(zip(fx.tolist(), fy.tolist())):
        ax.annotate(str(i + 1), (x, y), textcoords="offset points", xytext=(5, 4), fontsize=7, color="dimgrey")
    ax.set_xlabel(_axis_label(xk, xlo, xhi))
    ax.set_ylabel(_axis_label(yk, ylo, yhi))
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140)
    plt.close(fig)


def plot_3d(front, outfile, title="Frente de Pareto 3D", raw=None, ekey="Z4", zkey="Z5"):
    # Eje económico = Z4 (VAN total, el objetivo que realmente se optimiza en el
    # ε-constraint). NO Z1 (VAN max-min de la peor planta): Z1 es degenerado y toma
    # el piso del big-M (−M_pi_global≈−2e9) en la solución vacía, lo que aplastaba
    # la normalización del eje y producía la "L" (línea vertical + horizontal).
    # Eje ambiental = Z2 (GEI). Tercer eje = zkey: por defecto Z5 (output), el
    # tercer objetivo del ε-constraint. La cobertura Z3 quedó como diagnóstico y ya
    # no es eje del frente; para inspeccionarla, llamar con zkey="Z3".
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    basis = raw if raw else front
    z1lo, z1hi, _ = _minmax([p[ekey] for p in basis])
    z2lo, z2hi, _ = _minmax([p["Z2"] for p in basis])
    z3lo, z3hi, _ = _minmax([p[zkey] for p in basis])

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    if raw:
        ax.scatter(_scale([p[ekey] for p in raw], z1lo, z1hi),
                   _scale([p["Z2"] for p in raw], z2lo, z2hi),
                   _scale([p[zkey] for p in raw], z3lo, z3hi),
                   c="lightgrey", s=16, depthshade=True)
    ax.scatter(_scale([p[ekey] for p in front], z1lo, z1hi),
               _scale([p["Z2"] for p in front], z2lo, z2hi),
               _scale([p[zkey] for p in front], z3lo, z3hi),
               c="steelblue", s=50, edgecolors="k", linewidths=0.5, depthshade=True)
    ax.set_xlabel(_axis_label(ekey, z1lo, z1hi), fontsize=7)
    ax.set_ylabel(_axis_label("Z2", z2lo, z2hi), fontsize=7)
    ax.set_zlabel(_axis_label(zkey, z3lo, z3hi), fontsize=7)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140)
    plt.close(fig)


def plot_territory(front, geo, outfile):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    if not geo or not geo.get("coords_i"):
        print("  (datos geográficos no disponibles; omito plot_territory)")
        return

    coords_i = {int(k): v for k, v in geo["coords_i"].items()}
    coords_j = {int(k): v for k, v in geo["coords_j"].items()}
    coords_k = {int(k): v for k, v in geo["coords_k"].items()}
    A_JK = [(r[0], r[1]) for r in geo.get("A_JK", [])]

    freq = {}
    for p in front:
        for k in p["opened"]:
            freq[int(k)] = freq.get(int(k), 0) + 1
    max_freq = max(freq.values()) if freq else 1

    fig, ax = plt.subplots(figsize=(10, 9))
    xi = [v[0] for v in coords_i.values()]
    yi = [v[1] for v in coords_i.values()]
    wi = [v[2] for v in coords_i.values()]
    sc = ax.scatter(xi, yi, c=wi, cmap="YlOrRd", s=10, alpha=0.65, zorder=2)
    fig.colorbar(sc, ax=ax, label="Residuos anuales (kg/año)")

    open_ks = set(freq.keys())
    for (j, k) in A_JK:
        if k in open_ks and j in coords_j and k in coords_k:
            x0, y0 = coords_j[j]
            x1, y1 = coords_k[k][0], coords_k[k][1]
            ax.plot([x0, x1], [y0, y1], color="steelblue", linewidth=0.8, alpha=0.35, zorder=1)

    xj = [v[0] for v in coords_j.values()]
    yj = [v[1] for v in coords_j.values()]
    ax.scatter(xj, yj, marker="^", s=90, color="royalblue", edgecolors="k", linewidths=0.8, zorder=4)

    for k, v in coords_k.items():
        x, y = v[0], v[1]
        f = freq.get(k, 0)
        color = plt.cm.Greens(0.4 + 0.6 * f / max_freq) if f > 0 else "lightgrey"
        size = 60 + 140 * (f / max_freq) if f > 0 else 30
        ax.scatter(x, y, marker="s", s=size, color=color, edgecolors="k", linewidths=0.8, zorder=5)
        if f > 0:
            ax.annotate(f"k{k} ({f}/{len(front)})", (x, y), textcoords="offset points",
                        xytext=(5, 4), fontsize=6, color="darkgreen", zorder=6)

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="lightyellow", markeredgecolor="k", markersize=8),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="darkred", markeredgecolor="k", markersize=8),
        Line2D([0], [0], marker="^", color="w", markerfacecolor="royalblue", markeredgecolor="k", markersize=9),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="limegreen", markeredgecolor="k", markersize=9),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="lightgrey", markeredgecolor="k", markersize=9),
    ]
    ax.set_xlabel("X (L93)")
    ax.set_ylabel("Y (L93)")
    ax.set_title(f"Mapa territorial  |  D_BAR_JK={D_BAR_JK} km  |  {len(front)} soluciones Pareto")
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140)
    plt.close(fig)


def make_plots(front, outdir, raw=None, geo=None):
    if not front:
        print("  (sin puntos Pareto: nada que graficar)")
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        print("  matplotlib no disponible; omito figuras.")
        return
    except Exception:
        pass
    os.makedirs(outdir, exist_ok=True)
    # Excluye la solución trivial vacía (n_open=0): tiene Z1=−big-M y Z4=0 (no
    # invertir), y ensucia la normalización de todos los ejes económicos.
    fp = [p for p in front if p.get("n_open", 0) > 0]
    rp = [p for p in (raw or []) if p.get("n_open", 0) > 0] if raw else None
    if not fp:
        print("  (solo la solución vacía en el frente: nada informativo que graficar)")
        fp = front
    # Eje económico = Z4 (VAN total, objetivo real del ε-constraint). Z1 (max-min)
    # es degenerado y producía la "L"; se conserva en el JSON/CSV para diagnóstico.
    plot_2d(fp, "Z4", "Z2", os.path.join(outdir, "pareto_Z4_Z2.png"),
            "Pareto: económico (Z4=VAN total) vs ambiental (Z2)", raw=rp)
    plot_2d(fp, "Z3", "Z4", os.path.join(outdir, "pareto_Z3_Z4.png"),
            "Pareto: cobertura (Z3) vs económico (Z4=VAN total)", raw=rp)
    plot_2d(fp, "Z2", "Z3", os.path.join(outdir, "pareto_Z2_Z3.png"),
            "Pareto: ambiental (Z2) vs cobertura (Z3)", raw=rp)
    # Nuevos pares con Z5 (output)
    plot_2d(fp, "Z4", "Z5", os.path.join(outdir, "pareto_Z4_Z5.png"),
            "Pareto: económico (Z4=VAN total) vs output (Z5)", raw=rp)
    plot_2d(fp, "Z2", "Z5", os.path.join(outdir, "pareto_Z2_Z5.png"),
            "Pareto: ambiental (Z2) vs output (Z5)", raw=rp)
    plot_2d(fp, "Z3", "Z5", os.path.join(outdir, "pareto_Z3_Z5.png"),
            "Pareto: cobertura (Z3) vs output (Z5)", raw=rp)
    plot_3d(fp, os.path.join(outdir, "pareto_3d.png"),
            title="Frente de Pareto 3D (VAN · GEI · output)", raw=rp, ekey="Z4", zkey="Z5")
    plot_territory(fp, geo or {}, os.path.join(outdir, "territory.png"))
    print(f"  figuras → {outdir}/  (frente=Z4·Z2·Z5; Z3=cobertura solo diagnóstico)")


# ===========================================================================
# MAIN
# ===========================================================================
def _rng_from_payoff(payoff, constrained):
    """Reconstruye los rangos [min,max] por objetivo a partir de la tabla de pagos
    guardada en econ_results.json (payload['payoff']) — misma fórmula que
    lexicographic_payoff_table. Evita re-correr el payoff en el modo --time-curve.

    La clave 'payoff' del JSON no cambió al pasar al payoff lexicográfico (mismo
    formato {primary: {obj: valor}}), sólo cambió CÓMO se llenan las celdas: ahora
    cada fila es el punto final de una cadena lexicográfica, no un óptimo aislado.
    Un econ_results.json viejo (payoff mono-objetivo) sigue siendo legible, pero sus
    rangos son los inflados de antes — re-correr la grilla si se quiere consistencia."""
    rng = {}
    for k in constrained:
        vals = [c[k] for c in payoff.values() if isinstance(c, dict) and k in c]
        if vals:
            rng[k] = (float(min(vals)), float(max(vals)))
    return rng


def main():
    import argparse
    parser = argparse.ArgumentParser(description=f"sensibilidad.py — escenario: {SCENARIO}")
    parser.add_argument("--plot", action="store_true", help="Solo graficar desde resultados guardados")
    parser.add_argument("--time-curve", dest="time_curve", action="store_true",
                        help="Curva tiempo-vs-gap sobre las celdas problemáticas del último run")
    parser.add_argument("--results", default=None, metavar="PATH", help="Ruta al JSON de resultados")
    parser.add_argument("--outdir", default=None, metavar="DIR", help="Carpeta de salida")
    parser.add_argument("--zone", default=ZONE, metavar="Z", help="Zona a resolver")
    parser.add_argument("--list", action="store_true", help="Lista escenarios y zonas disponibles")
    args = parser.parse_args()

    if args.list:
        print(f"\nEscenarios disponibles (activo: {SCENARIO}):\n")
        for name, cfg in _SCENARIOS.items():
            active = " <-- ACTIVO" if name == SCENARIO else ""
            print(f"  {name:8s}  D_IJ={cfg['D_BAR_IJ']}  D_JK={cfg['D_BAR_JK']}  lca={str(cfg['LCA_CREDIT']):5s}{active}")
        man = os.path.join(INSTANCES_DIR, "manifest.json")
        if os.path.exists(man):
            with open(man, encoding="utf-8") as fh:
                zinfo = json.load(fh)["zones"]
            print(f"\nZonas disponibles (activa: {ZONE}):\n")
            for z in zinfo:
                act = " <-- ACTIVA" if z["name"] == ZONE else ""
                print(f"  {z['name']:14s}  I={z['n_I']:5d}  J={z['n_J']:3d}  K={z['n_K']:3d}  w={z['w_total']/1e6:.2f}M kg/año{act}")
        else:
            print("\n  (sin caché; corre: python build_instances.py)")
        print()
        return

    zone = args.zone
    outdir = args.outdir or os.path.join(
        REPO_ROOT, "results", f"results_sens_{zone}_{_S['tag']}")
    print(f"\n=== Escenario: {SCENARIO} ({_S['tag']}) | zona: {zone} | outdir: {outdir} ===")

    if args.plot:
        print("== Modo visualización (sin solver) ==")
        front, raw, geo, _ = load_frontier(args.results or os.path.join(outdir, "econ_results.json"))
        make_plots(front, outdir=outdir, raw=raw, geo=geo)
        return

    if args.time_curve:
        print("== Modo curva de optimalidad (tiempo vs gap sobre celdas problemáticas) ==")
        results_path = args.results or os.path.join(outdir, "econ_results.json")
        if not os.path.exists(results_path):
            print(f"No existe {results_path}. Corre primero la grilla (sin --time-curve).")
            return
        with open(results_path, encoding="utf-8") as fh:
            data = json.load(fh)
        raw = data.get("raw", [])
        payoff = data.get("payoff", {})
        rng = _rng_from_payoff(payoff, ECON_CONSTRAINED_ALL)
        if not rng:
            print("No se pudieron reconstruir los rangos del payoff desde el JSON. Aborto.")
            return
        # Un econ_results.json ANTERIOR al payoff lexicográfico puede no tener fila
        # para todos los objetivos de ECON_CONSTRAINED_ALL (p.ej. sin Z1). Se barren
        # sólo los que el JSON puede sostener: sin este filtro, run_time_curve haría
        # rng[k] -> KeyError sobre un objetivo ausente de la tabla guardada.
        tc_constrained = [k for k in ECON_CONSTRAINED_ALL if k in rng]
        if len(tc_constrained) < len(ECON_CONSTRAINED_ALL):
            falta = [k for k in ECON_CONSTRAINED_ALL if k not in rng]
            print(f"  ⚠ payoff guardado sin fila(s) para {falta} (resultado previo al payoff "
                  f"lexicográfico): la curva se corre sobre {tc_constrained}")
        prob = summarize_gaps(raw, tc_constrained)
        if not prob:
            print("\nSin celdas problemáticas: nada que barrer. Fin.")
            return
        print(f"\nCargando zona '{zone}' desde caché ({INSTANCES_DIR}) ...")
        D = load_instance(zone)
        if D is None:
            print(f"No se pudo cargar la instancia para la zona '{zone}'. Abortando.")
            return
        m = build_model(D)
        solver = make_solver()
        run_time_curve(m, solver, rng, prob, outdir, constrained=tc_constrained)
        return

    print(f"Cargando zona '{zone}' desde caché ({INSTANCES_DIR}) ...")
    D = load_instance(zone)
    if D is None:
        print(f"No se pudo cargar la instancia para la zona '{zone}'. Abortando.")
        return
    print(f"  |I|={len(D['I'])}  |J|={len(D['J'])}  |K|={len(D['K'])}  |M|={len(D['M'])}  |T|={len(D['TT'])}  "
          f"|A_IJ|={len(D['A_IJ'])}  |A_JK|={len(D['A_JK'])}  beta_min={D['beta_min']}")

    t0 = time.time()
    m = build_model(D)
    nv = sum(1 for _ in m.component_data_objects(pyo.Var, active=True))
    nc = sum(1 for _ in m.component_data_objects(pyo.Constraint, active=True))
    n_bin = len(D["K"]) * len(D["TT"])
    print(f"  modelo en {time.time()-t0:.1f}s: {nv:,} vars ({n_bin:,} binarias), {nc:,} restricciones")

    solver = make_solver()
    print("\n== Payoff table (lexicográfica, Mavrotas 2009 Sec. 3.1) ==")
    # Grilla pruebafull222: la tabla de pagos se construye SÓLO sobre los objetivos
    # que entran al ε-constraint — Z4 (VAN total, primario) + Z2 (GEI) y Z3
    # (cobertura). Z1 (VAN max-min) y Z5 (output) quedan COMENTADOS como objetivos:
    # su rango era degenerado/alineado y colapsaba el frente. Se siguen calculando y
    # reportando en CADA punto como columnas diagnósticas (econ_results.json,
    # soluciones.csv, plots), pero no consumen filas de payoff ni ejes de grilla.
    cells, rng, anchors = lexicographic_payoff_table(m, solver, objs=("Z4", "Z2", "Z3"))
    # cells, rng, anchors = lexicographic_payoff_table(m, solver, objs=("Z1", "Z2", "Z3", "Z4", "Z5"))  # 5-obj completo

    if RUN_ECONSTRAINT:
        front, raw, cells = epsilon_constraint(m, solver, cells, rng, seed_pts=anchors,
                                               constrained=ECON_CONSTRAINED_ALL,
                                               npts_map=ECON_NPTS_MAP)
        geo = _build_geo(D)
        save_frontier(front, raw, cells, outdir=outdir, geo=geo)
        save_detailed_csv(front, D, outdir)
        # Auditoría de gap: lista las celdas que no cerraron (entrada de --time-curve).
        summarize_gaps(raw)
        if MAKE_PLOTS:
            make_plots(front, outdir=outdir, raw=raw, geo=geo)


if __name__ == "__main__":
    main()
