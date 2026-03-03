param (
    [string]$Config = ".\configs\datasets_phase1b_densegrid_cleanradio.yaml",
    [string]$OutDir = ".\results\universality_run"
)

Write-Host "=== Protocol D – Universality Runner ==="

# Set PYTHONPATH
$env:PYTHONPATH = (Resolve-Path .\src).Path

# Optional raw sanity check (if script exists)
if (Test-Path ".\scripts\sanity_check_raw.py") {
    python .\scripts\sanity_check_raw.py --config $Config
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Raw file check failed. Aborting."
        exit 1
    }
}

# Run pipeline
python -m protocol_d.pipeline_universality `
    --config $Config `
    --out $OutDir

Write-Host "=== Run completed ==="