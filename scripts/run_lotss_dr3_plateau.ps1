# scripts/run_lotss_dr3_plateau.ps1
# LoTSS DR3 Test 1 (plateau) runner for Windows PowerShell
# - Uses Windows Python (same interpreter as 'python -m pip')
# - Avoids Git Bash / WSL conflicts
# - Assumes FITS is already downloaded, but can download if missing

$ErrorActionPreference = "Stop"

function Info($msg) { Write-Host "[INFO] $msg" }
function Warn($msg) { Write-Host "[WARN] $msg" -ForegroundColor Yellow }
function Fail($msg) { Write-Host "[ERROR] $msg" -ForegroundColor Red; exit 1 }

# Resolve repo root from this script location
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path

$FitsPath = Join-Path $RepoRoot "data\raw\LoTSS_DR3_v1.0.srl.fits"
$OutPath  = Join-Path $RepoRoot "results\lotss_dr3\plateau_dr3_counts_fit64.csv"
$ScriptPy = Join-Path $RepoRoot "src\lotss_dr3_plateau.py"

New-Item -ItemType Directory -Force -Path (Split-Path $FitsPath) | Out-Null
New-Item -ItemType Directory -Force -Path (Split-Path $OutPath)  | Out-Null

# Check Python availability
$pythonCmd = Get-Command python -ErrorAction SilentlyContinue
if (-not $pythonCmd) {
  Fail "python not found on PATH. Install Python and ensure 'python' works in PowerShell."
}

Info ("Python: " + (& python -c "import sys; print(sys.version.split()[0]); print(sys.executable)" | Out-String).Trim())

# Ensure deps (use same python)
$depsOk = $true
try {
  & python -c "import astropy; from astropy_healpix import healpy as hp; print('deps_ok')" | Out-Null
} catch {
  $depsOk = $false
}

if (-not $depsOk) {
  Warn "Missing deps (astropy / astropy-healpix). Installing into current user site-packages..."
  & python -m pip install --user --upgrade pip setuptools wheel
  & python -m pip install --user astropy astropy-healpix
}

# Download FITS if missing (optional)
if (-not (Test-Path $FitsPath)) {
  $url = "https://lofar-surveys.org/public/DR3/catalogues/LoTSS_DR3_v1.0.srl.fits"
  Info "FITS not found. Downloading from: $url"
  try { [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12 } catch {}
  try {
    Invoke-WebRequest -Uri $url -OutFile $FitsPath
  } catch {
    Fail "Download failed. Please download manually to: $FitsPath"
  }
} else {
  Info "Using existing FITS: $FitsPath"
}

# Run (absolute PYTHONPATH so it works from any working directory)
$env:PYTHONPATH = (Join-Path $RepoRoot "src")

Info "Running LoTSS DR3 plateau (fit_nside=64). Output -> $OutPath"

& python $ScriptPy `
  --fits $FitsPath `
  --out  $OutPath `
  --nside 128 `
  --fit_nside 64 `
  --flux_cuts_mjy "10,20,30,40,60,80" `
  --bmins_deg "10,20,30" `
  --footprint_flux_mjy 80

Info "Done."
Info "Output: $OutPath"
