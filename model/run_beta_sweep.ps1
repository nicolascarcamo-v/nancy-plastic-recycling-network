# Barrido de BETA_MIN (0.00 base, 0.05, 0.10, 0.20) sobre:
#   1. las 3 soluciones mono-objetivo  -> results/sensitivity_beta_mono/
#   2. la solucion de compromiso (21)  -> results/sensitivity_beta_compromise/
# y al final genera la figura figures/fig9_monthly_collection_beta.pdf.
#
# Uso (desde la raiz del repositorio):
#   powershell -ExecutionPolicy Bypass -File model/run_beta_sweep.ps1
# Reanuda solo: si se corta, volver a lanzarlo salta lo ya resuelto.

$ErrorActionPreference = "Continue"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$log = Join-Path $root "results/run_beta_sweep.log"

"=== Barrido beta iniciado: $(Get-Date) ===" | Tee-Object -FilePath $log -Append

"--- [1/3] mono-objetivo (sensibilidad_articulo.py --only beta) ---" | Tee-Object -FilePath $log -Append
python model/sensibilidad_articulo.py --only beta --outdir results/sensitivity_beta_mono --time-limit 600 2>&1 |
    Tee-Object -FilePath $log -Append

"--- [2/3] compromiso punto 21 (sensibilidad_compromiso.py --only beta) ---" | Tee-Object -FilePath $log -Append
python model/sensibilidad_compromiso.py --only beta --outdir results/sensitivity_beta_compromise --time-limit 1000 2>&1 |
    Tee-Object -FilePath $log -Append

"--- [3/3] figura fig9_monthly_collection_beta.pdf ---" | Tee-Object -FilePath $log -Append
python figures/scripts/fig9_monthly_collection_beta.py 2>&1 | Tee-Object -FilePath $log -Append

"=== Terminado: $(Get-Date) ===" | Tee-Object -FilePath $log -Append
