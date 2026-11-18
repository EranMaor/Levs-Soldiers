# Paper-scale block-error-rate run (Figs. 13-14). This takes days, not one night.
# Output is also written to results\full_bler.log.
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force -Path results | Out-Null
Write-Host "Full BLER started $(Get-Date -Format 'yyyy-MM-dd HH:mm'). Log: results\full_bler.log"
python -u plots.py bler --mode full 2>&1 | Tee-Object -FilePath results\full_bler.log
Write-Host "Full BLER finished $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
